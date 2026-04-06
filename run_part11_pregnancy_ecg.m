%% PART 1.1 - Abdominal ECG during pregnancy
% This script answers all items of Section 1.1.

clear; close all; clc;

[Y, fs, meta] = load_pregnancy_ecg();
leadNames = meta.leadNames;

%% 1) Plot the 8 leads
plot_multilead_ecg(Y, fs, leadNames, 'Pregnancy ECG - all 8 leads', size(Y,2)/fs);

disp('Observation: lead Abd4 is expected to be visibly noisier than the other abdominal leads.');

%% 2) Welch PSD of the noisy lead (lead 4)
noisyLead = Y(4,:);

[f4, P4, P4dB] = welch_psd_db(noisyLead, fs);
figure('Name','Lead Abd4 PSD','Color','w');
plot(f4, P4dB, 'k', 'LineWidth', 1.1); grid on;
xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('Welch PSD of noisy abdominal lead Abd4');

%% 3) Identify noise / interference intervals from the PSD
disp('Suggested interpretation from the PSD of Abd4:');
disp('- Very low-frequency content below about 0.5 Hz -> baseline wander');
disp('- Narrow interference around 50 Hz -> power-line contamination');
disp('- Broad high-frequency noise above about 40 Hz');

%% 4) Design filters
Hd = design_pregnancy_filters(fs);

lead4_filt = filtfilt(Hd.bHP, Hd.aHP, noisyLead);
lead4_filt = filtfilt(Hd.bBS, Hd.aBS, lead4_filt);
lead4_filt = filtfilt(Hd.bLP, Hd.aLP, lead4_filt);

%% 5) Apply filters and confirm effectiveness
[f4o, P4o, P4odB] = welch_psd_db(noisyLead, fs);
[f4f, P4f, P4fdB] = welch_psd_db(lead4_filt, fs);

t = (0:length(noisyLead)-1)/fs;
figure('Name','Lead Abd4 before/after filtering','Color','w');
subplot(2,1,1);
plot(t, noisyLead, 'Color', [0.3 0.3 0.3]); grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Noisy lead Abd4 - time domain');
subplot(2,1,2);
plot(t, lead4_filt, 'k'); grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Filtered lead Abd4 - time domain');

figure('Name','Lead Abd4 PSD before/after filtering','Color','w');
plot(f4o, P4odB, 'Color', [0.6 0.6 0.6], 'DisplayName', 'Original'); hold on;
plot(f4f, P4fdB, 'k', 'LineWidth', 1.1, 'DisplayName', 'Filtered');
grid on; xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('Lead Abd4 - PSD before and after filtering');
legend('Location','best');

results = struct();
results.lead4Original = noisyLead(:);
results.lead4Filtered = lead4_filt(:);
results.fs = fs;
results.fOriginal = f4o;
results.POriginal = P4o;
results.fFiltered = f4f;
results.PFiltered = P4f;

save('results_part11_pregnancy_ecg.mat', 'results');
disp('Saved: results_part11_pregnancy_ecg.mat');
