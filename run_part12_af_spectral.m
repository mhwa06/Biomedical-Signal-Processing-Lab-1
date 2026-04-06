%% PART 1.2 - Atrial fibrillation spectral analysis
% This script answers all items of Section 1.2 using the uploaded AF dataset.

clear; close all; clc;

[Xva, Xa, y, meta] = load_af_dataset();
fs = meta.fs;
v1 = meta.v1Idx;
pre = meta.precordialIdx;

[~, ~, nPatients] = size(Xva);

%% 1) DF and SC are computed by the helper function compute_df_sc.m

%% 2) Compute DF and SC in lead V1 of all ECG records
DF_Xva_V1 = NaN(nPatients,1);
SC_Xva_V1 = NaN(nPatients,1);
DF_Xa_V1  = NaN(nPatients,1);
SC_Xa_V1  = NaN(nPatients,1);

for p = 1:nPatients
    [DF_Xva_V1(p), SC_Xva_V1(p)] = compute_df_sc(Xva(:,v1,p), fs);
    [DF_Xa_V1(p),  SC_Xa_V1(p)]  = compute_df_sc(Xa(:, v1,p), fs);
end

%% 3) Boxplots: recorded data Xva versus estimated atrial activity Xa
figure('Name','V1 DF - Xva vs Xa','Color','w');
boxplot([DF_Xva_V1 DF_Xa_V1], 'Labels', {'Xva','Xa'});
ylabel('DF (Hz)'); grid on;
title('Dominant frequency in lead V1');

figure('Name','V1 SC - Xva vs Xa','Color','w');
boxplot([SC_Xva_V1 SC_Xa_V1], 'Labels', {'Xva','Xa'});
ylabel('SC (%)'); grid on;
title('Spectral concentration in lead V1');

%% 4) Classifier using DF and SC as features
res_V1_Xva = run_5fold_lda([DF_Xva_V1 SC_Xva_V1], y, '5-fold LDA using V1 features from Xva');
res_V1_Xa  = run_5fold_lda([DF_Xa_V1  SC_Xa_V1 ], y, '5-fold LDA using V1 features from Xa');

%% 5) Repeat parts 2 and 4 using averages over precordial leads V1:V6
nPre = numel(pre);
DF_Xva_pre = NaN(nPatients,nPre); SC_Xva_pre = NaN(nPatients,nPre);
DF_Xa_pre  = NaN(nPatients,nPre); SC_Xa_pre  = NaN(nPatients,nPre);

for p = 1:nPatients
    for k = 1:nPre
        ell = pre(k);
        [DF_Xva_pre(p,k), SC_Xva_pre(p,k)] = compute_df_sc(Xva(:,ell,p), fs);
        [DF_Xa_pre(p,k),  SC_Xa_pre(p,k)]  = compute_df_sc(Xa(:, ell,p), fs);
    end
end

DF_Xva_avg = mean(DF_Xva_pre, 2, 'omitnan');
SC_Xva_avg = mean(SC_Xva_pre, 2, 'omitnan');
DF_Xa_avg  = mean(DF_Xa_pre,  2, 'omitnan');
SC_Xa_avg  = mean(SC_Xa_pre,  2, 'omitnan');

res_AVG_Xva = run_5fold_lda([DF_Xva_avg SC_Xva_avg], y, '5-fold LDA using average V1:V6 features from Xva');
res_AVG_Xa  = run_5fold_lda([DF_Xa_avg  SC_Xa_avg ], y, '5-fold LDA using average V1:V6 features from Xa');

%% 6) Console conclusion support
fprintf('\nMedian values over labeled patients:\n');
idxLab = ~isnan(y);
fprintf('V1 Xva : median DF = %.3f Hz, median SC = %.3f %%\n', median(DF_Xva_V1(idxLab),'omitnan'), median(SC_Xva_V1(idxLab),'omitnan'));
fprintf('V1 Xa  : median DF = %.3f Hz, median SC = %.3f %%\n', median(DF_Xa_V1(idxLab),'omitnan'),  median(SC_Xa_V1(idxLab),'omitnan'));
fprintf('AVG Xva: median DF = %.3f Hz, median SC = %.3f %%\n', median(DF_Xva_avg(idxLab),'omitnan'), median(SC_Xva_avg(idxLab),'omitnan'));
fprintf('AVG Xa : median DF = %.3f Hz, median SC = %.3f %%\n', median(DF_Xa_avg(idxLab),'omitnan'),  median(SC_Xa_avg(idxLab),'omitnan'));

disp('Conclusion guide:');
disp('- Xa should usually have higher SC than Xva because the provided atrial estimate is less contaminated by QRST.');
disp('- Compare V1-only classification with averaged V1:V6 classification.');
disp('- State whether DF and SC alone give strong or only modest prediction of recurrence in this dataset.');

results = struct();
results.fs = fs;
results.labels = y;
results.DF_Xva_V1 = DF_Xva_V1;
results.SC_Xva_V1 = SC_Xva_V1;
results.DF_Xa_V1  = DF_Xa_V1;
results.SC_Xa_V1  = SC_Xa_V1;
results.DF_Xva_avg = DF_Xva_avg;
results.SC_Xva_avg = SC_Xva_avg;
results.DF_Xa_avg  = DF_Xa_avg;
results.SC_Xa_avg  = SC_Xa_avg;
results.res_V1_Xva = res_V1_Xva;
results.res_V1_Xa  = res_V1_Xa;
results.res_AVG_Xva = res_AVG_Xva;
results.res_AVG_Xa  = res_AVG_Xa;

save('results_part12_af_spectral.mat', 'results');
disp('Saved: results_part12_af_spectral.mat');
