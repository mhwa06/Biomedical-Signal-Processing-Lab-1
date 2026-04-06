function [idxAA, allDF, allSC] = select_aa_component(S, fs)
%SELECT_AA_COMPONENT Select the atrial component among estimated sources.
%
% Criterion:
%   among components whose DF is in [3, 9] Hz, keep the one with highest SC.

    cfg = assignment_config();

    nComp = size(S,1);
    allDF = NaN(nComp,1);
    allSC = NaN(nComp,1);

    for k = 1:nComp
        [allDF(k), allSC(k)] = compute_df_sc(S(k,:), fs);
    end

    cand = find(allDF >= cfg.dfBand(1) & allDF <= cfg.dfBand(2));

    if isempty(cand)
        [~, idxAA] = max(allSC);
    else
        [~, iBest] = max(allSC(cand));
        idxAA = cand(iBest);
    end
end
