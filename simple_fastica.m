function [S, H, W] = simple_fastica(X, maxIter, tol)
%SIMPLE_FASTICA Small self-contained FastICA implementation.
%
% Inputs
%   X       : observations, sensors x samples
%   maxIter : maximum iterations
%   tol     : stopping tolerance
%
% Outputs
%   S       : estimated sources, unit power, sources x samples
%   H       : estimated mixing matrix so that X ~= H*S
%   W       : demixing matrix on whitened data

    if nargin < 2 || isempty(maxIter), maxIter = 1000; end
    if nargin < 3 || isempty(tol), tol = 1e-6; end

    [p, T] = size(X);
    Xc = X - mean(X,2);

    C = (Xc * Xc.') / T;
    [E, D] = eig(C);
    [d, ord] = sort(diag(D), 'descend');
    E = E(:, ord);
    d = d(:);

    keep = d > max(d) * 1e-12;
    E = E(:, keep);
    d = d(keep);

    whitening = diag(1 ./ sqrt(d + eps)) * E.';
    Z = whitening * Xc;

    nComp = size(Z,1);
    W = zeros(nComp);

    rng(1);

    for k = 1:nComp
        w = randn(nComp,1);
        w = w / norm(w);

        for it = 1:maxIter
            wOld = w;

            proj = w.' * Z;           % 1 x T
            g = tanh(proj);
            gp = 1 - g.^2;

            w = (Z * g.') / T - mean(gp) * w;

            if k > 1
                w = w - W(:,1:k-1) * (W(:,1:k-1).' * w);
            end

            w = w / (norm(w) + eps);

            if abs(abs(w.' * wOld) - 1) < tol
                break;
            end
        end

        W(:,k) = w;
    end

    B = W.' * whitening;      % demixing matrix on centered data
    Sraw = B * Xc;            % nComp x T

    pow = sqrt(mean(Sraw.^2, 2) + eps);
    S = Sraw ./ pow;

    Hraw = pinv(B);
    H = Hraw * diag(pow);
end
