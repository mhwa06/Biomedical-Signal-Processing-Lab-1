function cfg = assignment_config()
%ASSIGNMENT_CONFIG Central configuration for the whole assignment.
%
% IMPORTANT:
%   The AF .mat files do not contain the sampling frequency explicitly.
%   The default below is set to 250 Hz because the records contain 15000
%   samples, which then corresponds to about 60 seconds, a common setup.
%   If your professor stated another value, change ONLY cfg.fsAF here.

    cfg = struct();

    % Sampling rates
    cfg.fsPreg = 500;
    cfg.fsAF   = 250;    % <-- change only this line if needed

    % Lead indexing
    cfg.afLeadNames = {'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'};
    cfg.afV1Idx = 7;
    cfg.afPrecordialIdx = 7:12;
    cfg.pregAbdominalIdx = 1:5;
    cfg.pregThoracicIdx  = 6:8;

    % Welch settings for AF spectral analysis
    cfg.dfBand = [3 9];              % Hz
    cfg.scBandFactor = [0.82 1.17];  % SC integration band = [0.82*DF, 1.17*DF]
    cfg.welch.windowLength = 4096;
    cfg.welch.overlapFactor = 0.50;
    cfg.welch.nfft = 8192;

    % Pregnancy ECG filtering
    cfg.preg.hpCutoff = 0.5;   % Hz, baseline wander removal
    cfg.preg.notch = [49 51];  % Hz, power line interference around 50 Hz
    cfg.preg.lpCutoff = 40;    % Hz, high-frequency noise removal

    % MRANC defaults for real data
    cfg.mranc.muPreg = 0.01;
    cfg.mranc.muAF   = 0.01;
    cfg.mranc.MReal  = 0;
    cfg.mranc.convSeconds = 2; % discard first 2 s when "after convergence" is required

    % Synthetic ANC defaults
    cfg.synth.N = 6000;
    cfg.synth.inputSIRdB = -5;
    cfg.synth.baseMu = 0.01;
    cfg.synth.baseM = 5;

    % ICA defaults
    cfg.ica.maxIter = 1000;
    cfg.ica.tol = 1e-6;
    cfg.ica.nonlinearity = 'tanh';
end
