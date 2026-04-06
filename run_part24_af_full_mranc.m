%% PART 2.4 - Atrial activity extraction by MRANC for the full AF database
% For each patient:
%   - primary = V1
%   - references = V2:V6
%   - evaluate DF and SC only

clear; close all; clc;

cfg = assignment_config();
[Xva, Xa, y, meta] = load_af_dataset();
fs = meta.fs;

nPatients = size(Xva,3);
nConv = round(cfg.mranc.convSeconds * fs);

DF_mranc = NaN(nPatients,1);
SC_mranc = NaN(nPatients,1);
DF_stc   = NaN(nPatients,1);
SC_stc   = NaN(nPatients,1);

for p = 1:nPatients
    Xpre = squeeze(Xva(:, meta.precordialIdx, p));
    d = Xpre(:,1);
    Xref = Xpre(:,2:6);

    [e, ~] = mranc_lms(d, Xref, cfg.mranc.muAF, 0);
    idxSS = (nConv+1):length(e);

    [DF_mranc(p), SC_mranc(p)] = compute_df_sc(e(idxSS), fs);
    [DF_stc(p),   SC_stc(p)]   = compute_df_sc(Xa(:,meta.v1Idx,p), fs);
end

figure('Name','MRANC vs STC - DF','Color','w');
boxplot([DF_mranc DF_stc], 'Labels', {'MRANC','STC'});
ylabel('DF (Hz)'); grid on;
title('Full database - dominant frequency');

figure('Name','MRANC vs STC - SC','Color','w');
boxplot([SC_mranc SC_stc], 'Labels', {'MRANC','STC'});
ylabel('SC (%)'); grid on;
title('Full database - spectral concentration');

fprintf('\nFull-database comparison over labeled patients:\n');
idxLab = ~isnan(y);
fprintf('MRANC: median DF = %.3f Hz, median SC = %.3f %%\n', median(DF_mranc(idxLab),'omitnan'), median(SC_mranc(idxLab),'omitnan'));
fprintf('STC  : median DF = %.3f Hz, median SC = %.3f %%\n', median(DF_stc(idxLab),'omitnan'),   median(SC_stc(idxLab),'omitnan'));

disp('Conclusion guide:');
disp('- If SC is higher, the atrial estimate is usually more organized and less contaminated by residual ventricular activity.');
disp('- Compare whether MRANC gives DF values comparable to STC and whether its SC is lower, similar, or higher.');

results = struct();
results.fs = fs;
results.labels = y;
results.DF_mranc = DF_mranc;
results.SC_mranc = SC_mranc;
results.DF_stc = DF_stc;
results.SC_stc = SC_stc;

save('results_part24_af_full_mranc.mat', 'results');
disp('Saved: results_part24_af_full_mranc.mat');

%% BONUS (optional): recurrence prediction from MRANC features
% Uncomment if you want the bonus.
% run_5fold_lda([DF_mranc SC_mranc], y, 'BONUS: 5-fold LDA using MRANC V1 features');
