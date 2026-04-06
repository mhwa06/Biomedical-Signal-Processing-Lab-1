function [f, Pxx, PxxdB] = welch_psd_db(x, fs)
%WELCH_PSD_DB Welch PSD using the assignment settings.

    cfg = assignment_config();

    x = x(:) - mean(x(:));
    N = length(x);

    wlen = min(cfg.welch.windowLength, N);
    wlen = max(wlen, min(N,256));

    win = hamming(wlen);
    noverlap = floor(cfg.welch.overlapFactor * wlen);
    nfft = max(cfg.welch.nfft, 2^nextpow2(wlen));

    [Pxx, f] = pwelch(x, win, noverlap, nfft, fs, 'onesided');
    PxxdB = 10*log10(Pxx + eps);
end
