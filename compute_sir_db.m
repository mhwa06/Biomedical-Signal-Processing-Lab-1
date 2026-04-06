function sir = compute_sir_db(sTrue, xObs)
%COMPUTE_SIR_DB Compute SIR in dB.
%
% sTrue is the desired signal.
% xObs is the observed or estimated signal = sTrue + interference/residual.

    sTrue = sTrue(:);
    xObs = xObs(:);
    err = xObs - sTrue;
    sir = 10*log10(var(sTrue) / (var(err) + eps));
end
