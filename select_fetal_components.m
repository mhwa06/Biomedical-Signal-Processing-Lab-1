function [idxFetal, scores, bpm, ratio] = select_fetal_components(H, S, fs)
%SELECT_FETAL_COMPONENTS Select likely fetal components from PCA/ICA outputs.
%
% Heuristic:
%   - fetal components should contribute more to abdominal than thoracic leads
%   - fetal heart rate is usually higher than the maternal one

    nComp = size(S,1);
    bpm = NaN(nComp,1);
    ratio = NaN(nComp,1);
    scores = NaN(nComp,1);

    for k = 1:nComp
        bpm(k) = estimate_heartbeat_bpm(S(k,:), fs);
        ratio(k) = (norm(H(1:5,k)) + eps) / (norm(H(6:8,k)) + eps);
        scores(k) = ratio(k) + 0.01 * bpm(k);
    end

    idx = find((ratio > 1.2) & (bpm > 100));

    if isempty(idx)
        [~, ord] = sort(scores, 'descend');
        idxFetal = ord(1:min(2,numel(ord)));
    else
        [~, ord] = sort(scores(idx), 'descend');
        idxFetal = idx(ord(1:min(2,numel(ord))));
    end
end
