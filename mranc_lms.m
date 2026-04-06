function [e, yhat, W] = mranc_lms(d, Xref, mu, M)
%MRANC_LMS Multi-reference adaptive noise canceling using LMS.
%
% Inputs
%   d     : primary input, N x 1
%   Xref  : reference matrix, N x R
%   mu    : LMS step size
%   M     : memory order
%
% Outputs
%   e     : canceller output (estimate of desired signal)
%   yhat  : estimated interference
%   W     : filter coefficients over time

    d = d(:);
    [N, R] = size(Xref);

    if size(Xref,1) ~= N
        error('Xref must have the same number of rows as d.');
    end

    L = R * (M + 1);
    e = zeros(N,1);
    yhat = zeros(N,1);
    W = zeros(L,N);

    w = zeros(L,1);

    for n = 1:N
        xvec = zeros(L,1);
        idx = 1;
        for r = 1:R
            for m = 0:M
                if n-m > 0
                    xvec(idx) = Xref(n-m, r);
                end
                idx = idx + 1;
            end
        end

        yhat(n) = w.' * xvec;
        e(n) = d(n) - yhat(n);
        w = w + 2 * mu * e(n) * xvec;
        W(:,n) = w;
    end
end
