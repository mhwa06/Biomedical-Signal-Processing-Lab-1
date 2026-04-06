%% PART 2.1 - Synthetic signals for MRANC
% NOTE:
% The exact synthetic setup of "Exercises 5 and 6 of the course" was not
% provided in the chat. This script gives a complete, runnable synthetic
% MRANC study that matches the questions: parameter effects and leakage.

clear; close all; clc;

cfg = assignment_config();
N = cfg.synth.N;
rng(1);

%% Build synthetic desired signal and correlated interferences
s = filter([1 0.5], [1 -1.4 0.65], randn(N,1));
s = s / std(s);

r1 = filter([1 0.3], [1 -0.6], randn(N,1)); r1 = r1/std(r1);
r2 = filter([1 -0.2 0.1], [1 -1.2 0.5], randn(N,1)); r2 = r2/std(r2);
r3 = filter([1 0.1], [1 -0.3], randn(N,1)); r3 = r3/std(r3);

h1 = [1 0.5 -0.2 0.1 0.05 0];
h2 = [0.8 -0.4 0.2 0.1 0 0];
h3 = [0.5 0.3 -0.2 0.1 -0.05 0];

v = conv(r1, h1, 'same') + conv(r2, h2, 'same') + conv(r3, h3, 'same');
v = v / std(v) * sqrt(10^(-cfg.synth.inputSIRdB/10));
d = s + v;

Xref = [r1 r2 r3];

figure('Name','Synthetic ANC signals','Color','w');
subplot(4,1,1); plot(d,'k'); grid on; title('Primary input d[n] = s[n] + v[n]');
subplot(4,1,2); plot(Xref(:,1),'k'); grid on; title('Reference 1');
subplot(4,1,3); plot(Xref(:,2),'k'); grid on; title('Reference 2');
subplot(4,1,4); plot(Xref(:,3),'k'); grid on; title('Reference 3'); xlabel('Samples');

%% 1) Effects of mu and M
tests = {
    0.005, 0;
    0.050, 0;
    0.010, 0;
    0.010, 5;
    0.010, 50
};

names = cell(size(tests,1),1);
SIRin = NaN(size(tests,1),1);
SIRout = NaN(size(tests,1),1);
SIRimp = NaN(size(tests,1),1);
nConv = NaN(size(tests,1),1);

for k = 1:size(tests,1)
    mu = tests{k,1};
    M = tests{k,2};

    [e, yhat] = mranc_lms(d, Xref, mu, M);
    nConv(k) = estimate_convergence_sample(s, e, 100);

    idxSS = nConv(k):N;
    SIRin(k) = compute_sir_db(s(idxSS), d(idxSS));
    SIRout(k) = compute_sir_db(s(idxSS), e(idxSS));
    SIRimp(k) = SIRout(k) - SIRin(k);
    names{k} = sprintf('mu=%.3f, M=%d', mu, M);

    if k <= 3
        figure('Name', names{k}, 'Color', 'w');
        subplot(3,1,1); plot(d,'Color',[0.6 0.6 0.6]); grid on; title(['Primary - ' names{k}]);
        subplot(3,1,2); plot(e,'k'); grid on; xline(nConv(k),'r--','Convergence');
        title('ANC output');
        subplot(3,1,3); plot(s,'b'); hold on; plot(e,'k'); grid on; xline(nConv(k),'r--');
        legend('True desired','ANC output'); title('Desired vs output');
    end
end

T = table(names, SIRin, SIRout, SIRimp, nConv, ...
    'VariableNames', {'Test','InputSIR_dB','OutputSIR_dB','SIRImprovement_dB','ConvergenceSample'});
disp(T);

disp('Conclusion guide for parameter effects:');
disp('- Larger mu usually speeds up convergence but may worsen steady-state performance or stability.');
disp('- Larger M can better model time-dispersive mixtures, but too large M increases the number of parameters and may slow convergence.');

%% 2) Leakage of desired signal into references (-10 dB)
XrefLeak = Xref + sqrt(0.1) * repmat(s, 1, size(Xref,2));

[eBase, ~] = mranc_lms(d, Xref, cfg.synth.baseMu, cfg.synth.baseM);
[eLeak, ~] = mranc_lms(d, XrefLeak, cfg.synth.baseMu, cfg.synth.baseM);

nConvBase = estimate_convergence_sample(s, eBase, 100);
nConvLeak = estimate_convergence_sample(s, eLeak, 100);

idxBase = nConvBase:N;
idxLeak = nConvLeak:N;

resLeak = struct();
resLeak.base.OutputSIR_dB = compute_sir_db(s(idxBase), eBase(idxBase));
resLeak.leak.OutputSIR_dB = compute_sir_db(s(idxLeak), eLeak(idxLeak));
resLeak.base.ConvergenceSample = nConvBase;
resLeak.leak.ConvergenceSample = nConvLeak;

figure('Name','Leakage effect in MRANC','Color','w');
subplot(2,1,1);
plot(eBase,'k'); hold on; plot(s,'b'); grid on; xline(nConvBase,'r--');
legend('ANC output - clean refs','True desired');
title('Without desired-signal leakage in references');
subplot(2,1,2);
plot(eLeak,'k'); hold on; plot(s,'b'); grid on; xline(nConvLeak,'r--');
legend('ANC output - contaminated refs','True desired');
title('With desired-signal leakage in references');

fprintf('\nLeakage experiment:\n');
fprintf('  Clean references     -> Output SIR = %.3f dB\n', resLeak.base.OutputSIR_dB);
fprintf('  Leaky references     -> Output SIR = %.3f dB\n', resLeak.leak.OutputSIR_dB);

disp('Conclusion guide for leakage:');
disp('- Leakage generally degrades ANC because the filter starts cancelling part of the desired signal.');
disp('- This produces signal distortion and usually lowers output SIR.');

results = struct();
results.table = T;
results.leakage = resLeak;
results.s = s;
results.d = d;
results.Xref = Xref;
results.XrefLeak = XrefLeak;

save('results_part21_synthetic_mranc.mat', 'results');
disp('Saved: results_part21_synthetic_mranc.mat');
