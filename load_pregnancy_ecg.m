function [Y, fs, meta] = load_pregnancy_ecg()
%LOAD_PREGNANCY_ECG Load the pregnant woman's ECG data.

    cfg = assignment_config();

    S = load('ecg_mother.mat');
    Y = S.Y;
    fs = cfg.fsPreg;

    meta = struct();
    meta.fs = fs;
    meta.abdominalIdx = cfg.pregAbdominalIdx;
    meta.thoracicIdx  = cfg.pregThoracicIdx;
    meta.leadNames = {'Abd1','Abd2','Abd3','Abd4','Abd5','Thor1','Thor2','Thor3'};

    fprintf('Loaded pregnancy ECG: %d leads x %d samples, fs = %g Hz\n', size(Y,1), size(Y,2), fs);
end
