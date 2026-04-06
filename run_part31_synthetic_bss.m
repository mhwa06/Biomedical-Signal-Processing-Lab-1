%% PART 3.1 - Synthetic signals for BSS
% NOTE:
% The exact "Exercises 5-7 of the course" were not provided. This script
% gives a full synthetic PCA/ICA study consistent with the questions.

clear; close all; clc;

cfg = assignment_config();
rng(2);
N = 6000;

%% Build 3 independent sources
s1 = filter([1 0.4], [1 -1.5 0.7], randn(N,1));  s1 = s1/std(s1);  % desired
t = (0:N-1).';
s2 = 2*mod(0.015*t,1)-1;                       s2 = s2/std(s2);
s3 = sign(sin(2*pi*0.03*(0:N-1))).';             s3 = s3/std(s3);

Strue = [s1.'; s2.'; s3.'];

%% 1) Reduce desired signal contribution in sensor 1 to 0.01
A1 = [0.01 0.7 0.5;
      0.8  1.0 0.2;
      0.5  0.3 1.0];
X1 = A1 * Strue;

[S_pca1, H_pca1] = bss_pca(X1);
[S_ica1, H_ica1] = simple_fastica(X1, cfg.ica.maxIter, cfg.ica.tol);

corr_pca = local_component_corr(S_pca1, s1);
corr_ica = local_component_corr(S_ica1, s1);
[bestCorrPCA, idxPCA] = max(corr_pca);
[bestCorrICA, idxICA] = max(corr_ica);

fprintf('\nExperiment 3.1-1 (desired contribution 0.01 in sensor 1)\n');
fprintf('  Best PCA correlation with desired source = %.3f (component %d)\n', bestCorrPCA, idxPCA);
fprintf('  Best ICA correlation with desired source = %.3f (component %d)\n', bestCorrICA, idxICA);

%% 2) Use observations with desired-signal leakage in reference-like sensors
A2 = [1.0  0.7 0.5;
      0.316 1.0 0.2;   % about -10 dB desired leakage
      0.316 0.3 1.0];
X2 = A2 * Strue;

[S_pca2, H_pca2] = bss_pca(X2);
[S_ica2, H_ica2] = simple_fastica(X2, cfg.ica.maxIter, cfg.ica.tol);

corr_pca2 = local_component_corr(S_pca2, s1);
corr_ica2 = local_component_corr(S_ica2, s1);
[bestCorrPCA2, idxPCA2] = max(corr_pca2);
[bestCorrICA2, idxICA2] = max(corr_ica2);

fprintf('\nExperiment 3.1-2 (leakage into reference-like sensors)\n');
fprintf('  Best PCA correlation with desired source = %.3f (component %d)\n', bestCorrPCA2, idxPCA2);
fprintf('  Best ICA correlation with desired source = %.3f (component %d)\n', bestCorrICA2, idxICA2);

%% Plots
figure('Name','Synthetic BSS - reduced desired contribution','Color','w');
subplot(3,1,1); plot(S_pca1(idxPCA,:),'k'); grid on; title('Best PCA component');
subplot(3,1,2); plot(S_ica1(idxICA,:),'k'); grid on; title('Best ICA component');
subplot(3,1,3); plot(s1,'b'); grid on; title('True desired source');

figure('Name','Synthetic BSS - leakage scenario','Color','w');
subplot(3,1,1); plot(S_pca2(idxPCA2,:),'k'); grid on; title('Best PCA component');
subplot(3,1,2); plot(S_ica2(idxICA2,:),'k'); grid on; title('Best ICA component');
subplot(3,1,3); plot(s1,'b'); grid on; title('True desired source');

disp('Conclusion guide:');
disp('- PCA uses only second-order information and may fail when the source of interest has weak variance contribution.');
disp('- ICA often separates non-Gaussian sources better than PCA.');
disp('- Leakage of the desired source into the other observations makes separation more difficult.');

results = struct();
results.bestCorrPCA_case1 = bestCorrPCA;
results.bestCorrICA_case1 = bestCorrICA;
results.bestCorrPCA_case2 = bestCorrPCA2;
results.bestCorrICA_case2 = bestCorrICA2;
save('results_part31_synthetic_bss.mat', 'results');
disp('Saved: results_part31_synthetic_bss.mat');


function c = local_component_corr(S, x)
    nComp = size(S,1);
    c = NaN(nComp,1);
    x = x(:);
    x = x - mean(x);
    for ii = 1:nComp
        s = S(ii,:).';
        s = s - mean(s);
        c(ii) = abs((s.'*x) / sqrt((s.'*s)*(x.'*x) + eps));
    end
end
