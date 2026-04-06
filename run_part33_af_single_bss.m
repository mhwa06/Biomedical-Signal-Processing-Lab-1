%% PART 3.3 - Atrial activity extraction by PCA and ICA for the first AF patient
% Observations = precordial leads V1:V6.

clear; close all; clc;

cfg = assignment_config();
[Xva, Xa, ~, meta] = load_af_dataset();
fs = meta.fs;

p = 1;
Xpre = squeeze(Xva(:, meta.precordialIdx, p));  % samples x 6
Xobs = Xpre.';                                  % sensors x samples

%% 1-2) Source separation via PCA and ICA
[S_pca, H_pca] = bss_pca(Xobs);
[S_ica, H_ica] = simple_fastica(Xobs, cfg.ica.maxIter, cfg.ica.tol);

%% 3a) Plot estimated sources
plot_multilead_ecg(S_pca, fs, arrayfun(@(k) sprintf('PCA %d',k), 1:size(S_pca,1), 'UniformOutput', false), ...
    'AF patient 1 - PCA sources', 8);

plot_multilead_ecg(S_ica, fs, arrayfun(@(k) sprintf('ICA %d',k), 1:size(S_ica,1), 'UniformOutput', false), ...
    'AF patient 1 - ICA sources', 8);

%% 3b-3d) Atrial-source identification
[idxAA_pca, DFpcaAll, SCpcaAll] = select_aa_component(S_pca, fs);
[idxAA_ica, DFicaAll, SCicaAll] = select_aa_component(S_ica, fs);

fprintf('\nAutomatic AA selection for patient 1:\n');
fprintf('  PCA selected component %d (DF = %.3f Hz, SC = %.3f %% )\n', idxAA_pca, DFpcaAll(idxAA_pca), SCpcaAll(idxAA_pca));
fprintf('  ICA selected component %d (DF = %.3f Hz, SC = %.3f %% )\n', idxAA_ica, DFicaAll(idxAA_ica), SCicaAll(idxAA_ica));

%% 3c) PSD of estimated sources
figure('Name','AF patient 1 - source PSDs (PCA)','Color','w');
for k = 1:size(S_pca,1)
    [f, P, PdB] = welch_psd_db(S_pca(k,:), fs);
    subplot(size(S_pca,1),1,k);
    plot(f, PdB, 'k'); grid on; xlim([0 25]);
    ylabel(sprintf('PCA %d',k));
    if k == 1, title('PCA source PSDs'); end
    if k == size(S_pca,1), xlabel('Frequency (Hz)'); end
end

figure('Name','AF patient 1 - source PSDs (ICA)','Color','w');
for k = 1:size(S_ica,1)
    [f, P, PdB] = welch_psd_db(S_ica(k,:), fs);
    subplot(size(S_ica,1),1,k);
    plot(f, PdB, 'k'); grid on; xlim([0 25]);
    ylabel(sprintf('ICA %d',k));
    if k == 1, title('ICA source PSDs'); end
    if k == size(S_ica,1), xlabel('Frequency (Hz)'); end
end

%% 3e) Contribution of the AA source to lead V1
aaV1_pca = H_pca(1,idxAA_pca) * S_pca(idxAA_pca,:);
aaV1_ica = H_ica(1,idxAA_ica) * S_ica(idxAA_ica,:);
v1Obs = Xpre(:,1).';
v1Stc = Xa(:, meta.v1Idx, p).';

%% 3f) Superimpose observed V1 and estimated AA contributions
t = (0:length(v1Obs)-1) / fs;
figure('Name','AF patient 1 - V1 and estimated AA contributions','Color','w');
subplot(3,1,1);
plot(t, v1Obs, 'Color',[0.6 0.6 0.6]); hold on; plot(t, aaV1_pca, 'k');
grid on; title('Observed V1 and PCA AA contribution'); legend('Observed V1','PCA AA');
subplot(3,1,2);
plot(t, v1Obs, 'Color',[0.6 0.6 0.6]); hold on; plot(t, aaV1_ica, 'k');
grid on; title('Observed V1 and ICA AA contribution'); legend('Observed V1','ICA AA');
subplot(3,1,3);
plot(t, v1Obs, 'Color',[0.6 0.6 0.6]); hold on; plot(t, v1Stc, 'k');
grid on; title('Observed V1 and STC atrial estimate'); legend('Observed V1','STC');
xlabel('Time (s)');

%% 3g) PSD comparison between observed V1 and estimated AA contributions
[f1, P1, P1dB] = welch_psd_db(v1Obs, fs);
[f2, P2, P2dB] = welch_psd_db(aaV1_pca, fs);
[f3, P3, P3dB] = welch_psd_db(aaV1_ica, fs);

figure('Name','AF patient 1 - PSD comparison at V1','Color','w');
plot(f1, P1dB, 'Color',[0.6 0.6 0.6], 'DisplayName','Observed V1'); hold on;
plot(f2, P2dB, 'k', 'DisplayName','PCA AA contribution');
plot(f3, P3dB, '--k', 'DisplayName','ICA AA contribution');
grid on; xlim([0 25]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('Observed V1 vs estimated atrial contribution at V1');
legend('Location','best');

%% 3h) DF and SC on the estimated atrial signal in lead V1
[DF_v1_pca, SC_v1_pca] = compute_df_sc(aaV1_pca, fs);
[DF_v1_ica, SC_v1_ica] = compute_df_sc(aaV1_ica, fs);
[DF_v1_stc, SC_v1_stc] = compute_df_sc(v1Stc, fs);

fprintf('\nSpectral indices on lead V1 atrial estimates (patient 1):\n');
fprintf('  PCA: DF = %.3f Hz, SC = %.3f %%\n', DF_v1_pca, SC_v1_pca);
fprintf('  ICA: DF = %.3f Hz, SC = %.3f %%\n', DF_v1_ica, SC_v1_ica);
fprintf('  STC: DF = %.3f Hz, SC = %.3f %%\n', DF_v1_stc, SC_v1_stc);

disp('Comparison guide with MRANC and STC:');
disp('- Compare which method gives the clearest narrow-band atrial peak in the AF band.');
disp('- Higher SC usually indicates a cleaner atrial estimate.');

results = struct();
results.patient = p;
results.fs = fs;
results.idxAA_pca = idxAA_pca;
results.idxAA_ica = idxAA_ica;
results.aaV1_pca = aaV1_pca;
results.aaV1_ica = aaV1_ica;
results.v1Stc = v1Stc;
results.DF_v1_pca = DF_v1_pca;
results.SC_v1_pca = SC_v1_pca;
results.DF_v1_ica = DF_v1_ica;
results.SC_v1_ica = SC_v1_ica;
results.DF_v1_stc = DF_v1_stc;
results.SC_v1_stc = SC_v1_stc;

save('results_part33_af_single_bss.mat', 'results');
disp('Saved: results_part33_af_single_bss.mat');
