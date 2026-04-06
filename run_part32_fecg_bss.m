%% PART 3.2 - Fetal ECG extraction by PCA and ICA
% Observations = the 8 pregnancy ECG channels.

clear; close all; clc;

cfg = assignment_config();
[Y, fs, meta] = load_pregnancy_ecg();

Xobs = Y;   % sensors x samples

%% 1) Perform source separation
[S_pca, H_pca] = bss_pca(Xobs);
[S_ica, H_ica] = simple_fastica(Xobs, cfg.ica.maxIter, cfg.ica.tol);

%% 2a) Plot sources
plot_multilead_ecg(S_pca, fs, arrayfun(@(k) sprintf('PCA %d',k), 1:size(S_pca,1), 'UniformOutput', false), ...
    'PCA sources - pregnancy ECG', 8);

plot_multilead_ecg(S_ica, fs, arrayfun(@(k) sprintf('ICA %d',k), 1:size(S_ica,1), 'UniformOutput', false), ...
    'ICA sources - pregnancy ECG', 8);

%% 2b) Visually identify fetal ECG sources (automatic aid)
[idxF_pca, scorePCA, bpmPCA, ratioPCA] = select_fetal_components(H_pca, S_pca, fs);
[idxF_ica, scoreICA, bpmICA, ratioICA] = select_fetal_components(H_ica, S_ica, fs);

fprintf('\nSuggested fetal components based on abdominal/thoracic ratio and heart rate:\n');
fprintf('  PCA candidates: %s\n', mat2str(idxF_pca(:).'));
fprintf('  ICA candidates: %s\n', mat2str(idxF_ica(:).'));

%% 2c) Compute fetal source contribution to all recordings
Xfetal_pca = H_pca(:, idxF_pca) * S_pca(idxF_pca, :);
Xfetal_ica = H_ica(:, idxF_ica) * S_ica(idxF_ica, :);

%% 2d) Plot abdominal signals and estimated fetal contributions
t = (0:size(Y,2)-1) / fs;
figure('Name','Pregnancy BSS - abdominal signals vs PCA fetal contribution','Color','w');
for k = 1:5
    subplot(5,1,k);
    plot(t, Y(k,:), 'Color',[0.6 0.6 0.6], 'DisplayName','Abdominal signal'); hold on;
    plot(t, Xfetal_pca(k,:), 'k', 'DisplayName','Estimated fetal contribution');
    grid on; ylabel(sprintf('Abd%d',k));
    if k == 1, title('PCA: abdominal signals vs estimated fetal contribution'); end
    if k == 5, xlabel('Time (s)'); end
    legend('Location','best');
end

figure('Name','Pregnancy BSS - abdominal signals vs ICA fetal contribution','Color','w');
for k = 1:5
    subplot(5,1,k);
    plot(t, Y(k,:), 'Color',[0.6 0.6 0.6], 'DisplayName','Abdominal signal'); hold on;
    plot(t, Xfetal_ica(k,:), 'k', 'DisplayName','Estimated fetal contribution');
    grid on; ylabel(sprintf('Abd%d',k));
    if k == 1, title('ICA: abdominal signals vs estimated fetal contribution'); end
    if k == 5, xlabel('Time (s)'); end
    legend('Location','best');
end

%% 2e) Compare spectra 0-50 Hz
figure('Name','Pregnancy BSS - PSD comparison PCA','Color','w');
for k = 1:5
    [f1, P1, P1dB] = welch_psd_db(Y(k,:), fs);
    [f2, P2, P2dB] = welch_psd_db(Xfetal_pca(k,:), fs);
    subplot(5,1,k);
    plot(f1, P1dB, 'Color',[0.6 0.6 0.6]); hold on;
    plot(f2, P2dB, 'k');
    grid on; xlim([0 50]);
    ylabel(sprintf('Abd%d',k));
    if k == 1, title('PCA: abdominal signals vs fetal contribution PSD'); end
    if k == 5, xlabel('Frequency (Hz)'); end
end

figure('Name','Pregnancy BSS - PSD comparison ICA','Color','w');
for k = 1:5
    [f1, P1, P1dB] = welch_psd_db(Y(k,:), fs);
    [f2, P2, P2dB] = welch_psd_db(Xfetal_ica(k,:), fs);
    subplot(5,1,k);
    plot(f1, P1dB, 'Color',[0.6 0.6 0.6]); hold on;
    plot(f2, P2dB, 'k');
    grid on; xlim([0 50]);
    ylabel(sprintf('Abd%d',k));
    if k == 1, title('ICA: abdominal signals vs fetal contribution PSD'); end
    if k == 5, xlabel('Frequency (Hz)'); end
end

%% 3) Compare with MRANC results
disp('Comparison guide with Part 2.2 MRANC:');
disp('- MRANC uses thoracic leads as references for maternal cancellation.');
disp('- PCA/ICA are blind separation methods using all 8 leads without explicit references.');
disp('- In your discussion, compare visual quality, fetal contribution clarity, and spectral content.');

results = struct();
results.fs = fs;
results.idxF_pca = idxF_pca;
results.idxF_ica = idxF_ica;
results.Xfetal_pca = Xfetal_pca;
results.Xfetal_ica = Xfetal_ica;
results.bpmPCA = bpmPCA;
results.bpmICA = bpmICA;
results.ratioPCA = ratioPCA;
results.ratioICA = ratioICA;
save('results_part32_fecg_bss.mat', 'results');

disp('Saved: results_part32_fecg_bss.mat');
