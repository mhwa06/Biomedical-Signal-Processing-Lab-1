%% PART 2.2 - Fetal ECG extraction by MRANC
% Abdominal leads are primaries; thoracic leads are references; M = 0.

clear; close all; clc;

cfg = assignment_config();
[Y, fs, meta] = load_pregnancy_ecg();

abd = cfg.pregAbdominalIdx;
thor = cfg.pregThoracicIdx;
mu = cfg.mranc.muPreg;
M = 0;
nConv = round(cfg.mranc.convSeconds * fs);

Xref = Y(thor,:).';
nAbd = numel(abd);

fetalEst = zeros(size(Y));

for k = 1:nAbd
    d = Y(abd(k), :).';
    [e, yhat] = mranc_lms(d, Xref, mu, M);
    fetalEst(abd(k), :) = e.';
end

%% (a) Plot inputs and outputs, zoom in first 2 seconds
t = (0:size(Y,2)-1) / fs;
tZoom = t <= 2;

figure('Name','FECG extraction by MRANC - first 2 seconds','Color','w');
for k = 1:nAbd
    subplot(nAbd,1,k);
    plot(t(tZoom), Y(abd(k),tZoom), 'Color',[0.6 0.6 0.6], 'DisplayName','Primary'); hold on;
    plot(t(tZoom), fetalEst(abd(k),tZoom), 'k', 'DisplayName','Output');
    grid on;
    ylabel(sprintf('Abd%d', k));
    if k == 1
        title('Primary abdominal leads vs MRANC outputs (first 2 s)');
    end
    if k == nAbd
        xlabel('Time (s)');
    end
    legend('Location','best');
end

%% (b) Compare PSDs between 0 and 50 Hz after convergence
figure('Name','FECG MRANC - PSD comparison','Color','w');
for k = 1:nAbd
    xIn = Y(abd(k), :);
    xOut = fetalEst(abd(k), nConv+1:end);

    [fIn, PIn, PIndB] = welch_psd_db(xIn, fs);
    [fOut, POut, POutdB] = welch_psd_db(xOut, fs);

    subplot(nAbd,1,k);
    plot(fIn, PIndB, 'Color',[0.6 0.6 0.6], 'DisplayName','Primary'); hold on;
    plot(fOut, POutdB, 'k', 'DisplayName','Output');
    grid on; xlim([0 50]);
    ylabel(sprintf('Abd%d', k));
    if k == 1
        title('Welch PSD before and after MRANC (after convergence)');
    end
    if k == nAbd
        xlabel('Frequency (Hz)');
    end
    legend('Location','best');
end

results = struct();
results.fs = fs;
results.Y = Y;
results.fetalEst = fetalEst;
results.nConv = nConv;
save('results_part22_fecg_mranc.mat', 'results');

disp('Saved: results_part22_fecg_mranc.mat');
disp('Interpretation guide: the MRANC output should suppress the large maternal component and make fetal beats more visible.');
