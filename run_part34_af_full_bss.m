%% PART 3.4 - Atrial activity extraction by PCA and ICA for the full AF database
% Automatic identification assumes a single atrial source.
% This script compares PCA, ICA, MRANC, and STC using DF and SC.

clear; close all; clc;

cfg = assignment_config();
[Xva, Xa, y, meta] = load_af_dataset();
fs = meta.fs;
nPatients = size(Xva,3);
nConv = round(cfg.mranc.convSeconds * fs);

DF_pca = NaN(nPatients,1); SC_pca = NaN(nPatients,1);
DF_ica = NaN(nPatients,1); SC_ica = NaN(nPatients,1);
DF_mranc = NaN(nPatients,1); SC_mranc = NaN(nPatients,1);
DF_stc = NaN(nPatients,1);  SC_stc = NaN(nPatients,1);

for p = 1:nPatients
    Xpre = squeeze(Xva(:, meta.precordialIdx, p));   % samples x 6
    Xobs = Xpre.';                                   % sensors x samples

    % PCA
    [S_pca, H_pca] = bss_pca(Xobs);
    [idxAA_pca, ~, ~] = select_aa_component(S_pca, fs);
    aaV1_pca = H_pca(1,idxAA_pca) * S_pca(idxAA_pca,:);
    [DF_pca(p), SC_pca(p)] = compute_df_sc(aaV1_pca, fs);

    % ICA
    [S_ica, H_ica] = simple_fastica(Xobs, cfg.ica.maxIter, cfg.ica.tol);
    [idxAA_ica, ~, ~] = select_aa_component(S_ica, fs);
    aaV1_ica = H_ica(1,idxAA_ica) * S_ica(idxAA_ica,:);
    [DF_ica(p), SC_ica(p)] = compute_df_sc(aaV1_ica, fs);

    % MRANC
    d = Xpre(:,1);
    Xref = Xpre(:,2:6);
    [e, ~] = mranc_lms(d, Xref, cfg.mranc.muAF, 0);
    idxSS = (nConv+1):length(e);
    [DF_mranc(p), SC_mranc(p)] = compute_df_sc(e(idxSS), fs);

    % STC (provided)
    [DF_stc(p), SC_stc(p)] = compute_df_sc(Xa(:,meta.v1Idx,p), fs);
end

figure('Name','Full database comparison - DF','Color','w');
boxplot([DF_pca DF_ica DF_mranc DF_stc], 'Labels', {'PCA','ICA','MRANC','STC'});
ylabel('DF (Hz)'); grid on;
title('Full database - dominant frequency');

figure('Name','Full database comparison - SC','Color','w');
boxplot([SC_pca SC_ica SC_mranc SC_stc], 'Labels', {'PCA','ICA','MRANC','STC'});
ylabel('SC (%)'); grid on;
title('Full database - spectral concentration');

idxLab = ~isnan(y);
fprintf('\nMedian DF over labeled patients:\n');
fprintf('  PCA   = %.3f Hz\n', median(DF_pca(idxLab),'omitnan'));
fprintf('  ICA   = %.3f Hz\n', median(DF_ica(idxLab),'omitnan'));
fprintf('  MRANC = %.3f Hz\n', median(DF_mranc(idxLab),'omitnan'));
fprintf('  STC   = %.3f Hz\n', median(DF_stc(idxLab),'omitnan'));

fprintf('\nMedian SC over labeled patients:\n');
fprintf('  PCA   = %.3f %%\n', median(SC_pca(idxLab),'omitnan'));
fprintf('  ICA   = %.3f %%\n', median(SC_ica(idxLab),'omitnan'));
fprintf('  MRANC = %.3f %%\n', median(SC_mranc(idxLab),'omitnan'));
fprintf('  STC   = %.3f %%\n', median(SC_stc(idxLab),'omitnan'));

disp('Conclusion guide:');
disp('- Compare organization/cleanliness through SC and physiological consistency through DF.');
disp('- In many AF applications, ICA tends to outperform PCA because independence is more appropriate than orthogonality for atrial/ventricular separation.');
disp('- Use your boxplots and median values to conclude which method is best on this dataset.');

results = struct();
results.fs = fs;
results.labels = y;
results.DF_pca = DF_pca;     results.SC_pca = SC_pca;
results.DF_ica = DF_ica;     results.SC_ica = SC_ica;
results.DF_mranc = DF_mranc; results.SC_mranc = SC_mranc;
results.DF_stc = DF_stc;     results.SC_stc = SC_stc;

save('results_part34_af_full_bss.mat', 'results');
disp('Saved: results_part34_af_full_bss.mat');

%% BONUS (optional): recurrence prediction
% Uncomment to test the bonus requested in the statement.
% run_5fold_lda([DF_pca SC_pca], y, 'BONUS: PCA features');
% run_5fold_lda([DF_ica SC_ica], y, 'BONUS: ICA features');
