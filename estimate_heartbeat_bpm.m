function bpm = estimate_heartbeat_bpm(x, fs)
%ESTIMATE_HEARTBEAT_BPM Estimate a heart rate from the dominant low-frequency peak.

    x = x(:) - mean(x(:));
    [f, Pxx] = welch_psd_db(x, fs);

    idx = (f >= 0.8) & (f <= 4.0);
    if ~any(idx)
        bpm = NaN;
        return;
    end

    [~, iMax] = max(Pxx(idx));
    fBand = f(idx);
    bpm = 60 * fBand(iMax);
end
