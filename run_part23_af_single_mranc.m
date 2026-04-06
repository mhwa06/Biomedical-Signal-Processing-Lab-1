%% PART 2.3 - Atrial activity extraction by MRANC for the first patient
% Primary input = V1, references = V2:V6, M = 0.

clear; close all; clc;

cfg = assignment_config();
[Xva, Xa, ~, meta] = load_af_dataset();
fs = meta.fs;

p = 1;
Xpre = squeeze(Xva(:, meta.precordialIdx, p));   % samples x 6
Xa_pre = squeeze(Xa(:, meta.precordialIdx, p));

%% 1) Plot the 6-lead precordial record
plot_multilead_ecg(Xpre, fs, meta.leadNames(meta.precordialIdx), 'Patient 1 - precordial leads V1:V6', 8);
disp('Observation: V1 is usually the lead where atrial activity is more visible.');

%% 2) MRANC: V1 primary, V2:V6 references
d = Xpre(:,1);
Xref = Xpre(:,2:6);
[e, yhat] = mranc_lms(d, Xref, cfg.mranc.muAF, 0);

nConv = round(cfg.mranc.convSeconds * fs);
idxSS = (nConv+1):length(e);

%% 3a) Time-domain comparison after convergence
t = (0:length(d)-1) / fs;
tWin = (t >= cfg.mranc.convSeconds) & (t <= cfg.mranc.convSeconds + 2);

figure('Name','AF single-patient MRANC - time domain','Color','w');
plot(t(tWin), d(tWin), 'Color',[0.6 0.6 0.6], 'DisplayName','Primary V1'); hold on;
plot(t(tWin), e(tWin), 'k', 'DisplayName','MRANC output');
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Patient 1 - V1 and MRANC output (2-second window after convergence)');
legend('Location','best');

%% 3b) PSD comparison 0-50 Hz after convergence
[fIn, PIn, PIndB] = welch_psd_db(d, fs);
[fOut, POut, POutdB] = welch_psd_db(e(idxSS), fs);

figure('Name','AF single-patient MRANC - PSD','Color','w');
plot(fIn, PIndB, 'Color',[0.6 0.6 0.6], 'DisplayName','Primary V1'); hold on;
plot(fOut, POutdB, 'k', 'DisplayName','MRANC output');
grid on; xlim([0 50]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('Patient 1 - PSD comparison');
legend('Location','best');

%% 3c) DF and SC on the output after convergence
[DF_mranc, SC_mranc] = compute_df_sc(e(idxSS), fs);

%% 3d) Compare with the provided STC signal in Xa
xaa_v1 = Xa_pre(:,1);
[DF_stc, SC_stc] = compute_df_sc(xaa_v1, fs);

fprintf('\nPatient 1 spectral parameters:\n');
fprintf('  MRANC output: DF = %.3f Hz, SC = %.3f %%\n', DF_mranc, SC_mranc);
fprintf('  STC (Xa V1):  DF = %.3f Hz, SC = %.3f %%\n', DF_stc, SC_stc);

results = struct();
results.patient = p;
results.fs = fs;
results.primaryV1 = d;
results.mrancOutput = e;
results.stcV1 = xaa_v1;
results.nConv = nConv;
results.DF_mranc = DF_mranc;
results.SC_mranc = SC_mranc;
results.DF_stc = DF_stc;
results.SC_stc = SC_stc;

save('results_part23_af_single_mranc.mat', 'results');
disp('Saved: results_part23_af_single_mranc.mat');
