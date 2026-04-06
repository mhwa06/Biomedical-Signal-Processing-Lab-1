function [Xva, Xa, y, meta] = load_af_dataset()
%LOAD_AF_DATASET Load the AF dataset from the uploaded MAT files.
%
% Outputs:
%   Xva  - recorded ECG, size = samples x leads x patients
%   Xa   - provided atrial estimate (STC), same size
%   y    - labels, length = patients
%          0 = non-recurrence, 1 = recurrence, NaN = unlabeled
%   meta - structure with lead indices and metadata

    cfg = assignment_config();

    A = load('Rva1.mat'); B = load('Rva2.mat'); C = load('Rva3.mat');
    D = load('Ra1.mat');  E = load('Ra2.mat');  F = load('Ra3.mat');
    G = load('indnonrecur.mat');
    H = load('indrecur.mat');

    Rva = cat(3, A.Rva1, B.Rva2, C.Rva3);   % leads x samples x patients
    Ra  = cat(3, D.Ra1,  E.Ra2,  F.Ra3);

    Xva = permute(Rva, [2 1 3]);            % samples x leads x patients
    Xa  = permute(Ra,  [2 1 3]);

    nPatients = size(Xva, 3);
    y = NaN(nPatients, 1);

    indnon = double(G.indnonrecur(:));
    indrec = double(H.indrecur(:));

    y(indnon) = 0;
    y(indrec) = 1;

    meta = struct();
    meta.fs = cfg.fsAF;
    meta.leadNames = cfg.afLeadNames;
    meta.v1Idx = cfg.afV1Idx;
    meta.precordialIdx = cfg.afPrecordialIdx;
    meta.patientIDs = (1:nPatients).';
    meta.labeledIdx = find(~isnan(y));

    fprintf('Loaded AF dataset: %d samples x %d leads x %d patients\n', size(Xva,1), size(Xva,2), size(Xva,3));
    fprintf('Labeled patients: %d / %d\n', numel(meta.labeledIdx), nPatients);
    fprintf('AF sampling frequency currently set to %g Hz (see assignment_config.m)\n', meta.fs);
end
