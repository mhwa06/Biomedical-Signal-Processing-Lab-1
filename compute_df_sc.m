function [DF, SC, f, Pxx] = compute_df_sc(x, fs)
%COMPUTE_DF_SC Compute Dominant Frequency and Spectral Concentration.
%
% DF = position of the dominant PSD peak in [3, 9] Hz.
% SC = 100 * power in [0.82*DF, 1.17*DF] / total one-sided power.

    cfg = assignment_config();

    x = x(:) - mean(x(:));
    [f, Pxx] = welch_psd_db(x, fs);
    Pxx = Pxx(:);
    f = f(:);

    band = (f >= cfg.dfBand(1)) & (f <= cfg.dfBand(2));
    if ~any(band)
        DF = NaN; SC = NaN; return;
    end

    [~, iPeak] = max(Pxx(band));
    idxBand = find(band);
    idxPeak = idxBand(iPeak);
    DF = f(idxPeak);

    scLow = cfg.scBandFactor(1) * DF;
    scHigh = cfg.scBandFactor(2) * DF;

    idxSC = (f >= scLow) & (f <= scHigh);
    totalPow = trapz(f, Pxx);
    scPow = trapz(f(idxSC), Pxx(idxSC));

    SC = 100 * scPow / (totalPow + eps);
end
