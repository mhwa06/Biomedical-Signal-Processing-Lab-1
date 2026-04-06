function nConv = estimate_convergence_sample(sTrue, sEst, fs)
%ESTIMATE_CONVERGENCE_SAMPLE Estimate a convergence sample for ANC output.
%
% The criterion is based on a running output SIR reaching within 1 dB
% of the final steady-state SIR.

    if nargin < 3
        fs = 1;
    end

    sTrue = sTrue(:);
    sEst = sEst(:);
    N = length(sTrue);

    win = max(round(0.5*fs), 100);
    if N < 3*win
        nConv = min(round(0.2*N), N);
        return;
    end

    outSIR = NaN(N-win+1,1);
    for n = 1:(N-win+1)
        segTrue = sTrue(n:n+win-1);
        segErr = sEst(n:n+win-1) - segTrue;
        outSIR(n) = 10*log10(var(segTrue) / (var(segErr) + eps));
    end

    finalSIR = mean(outSIR(max(1,end-10):end), 'omitnan');
    idx = find(outSIR >= finalSIR - 1, 1, 'first');

    if isempty(idx)
        nConv = min(round(0.2*N), N);
    else
        nConv = idx;
    end
end
