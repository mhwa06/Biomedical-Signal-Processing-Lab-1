function res = run_5fold_lda(X, y, tag)
%RUN_5FOLD_LDA 5-fold stratified cross-validation with a simple LDA classifier.
%
% Inputs
%   X   : feature matrix, N x D
%   y   : labels in {0,1,NaN}
%   tag : label used for printing
%
% Output
%   res : structure with predictions and metrics

    if nargin < 3
        tag = '5-fold LDA';
    end

    valid = ~isnan(y) & all(~isnan(X),2);
    X = X(valid,:);
    y = y(valid);

    rng(1);

    idx0 = find(y == 0);
    idx1 = find(y == 1);

    idx0 = idx0(randperm(numel(idx0)));
    idx1 = idx1(randperm(numel(idx1)));

    folds = cell(5,1);
    for k = 1:5
        folds{k} = [];
    end
    for i = 1:numel(idx0)
        folds{mod(i-1,5)+1} = [folds{mod(i-1,5)+1}; idx0(i)];
    end
    for i = 1:numel(idx1)
        folds{mod(i-1,5)+1} = [folds{mod(i-1,5)+1}; idx1(i)];
    end

    yhat = NaN(size(y));
    score = NaN(size(y));

    for k = 1:5
        testIdx = sort(folds{k});
        trainMask = true(size(y));
        trainMask(testIdx) = false;
        trainIdx = find(trainMask);

        Xtr = X(trainIdx,:);
        ytr = y(trainIdx);
        Xte = X(testIdx,:);

        mu = mean(Xtr, 1);
        sigma = std(Xtr, 0, 1);
        sigma(sigma < eps) = 1;

        Xtrz = (Xtr - mu) ./ sigma;
        Xtez = (Xte - mu) ./ sigma;

        [w, c] = local_train_lda(Xtrz, ytr);
        score(testIdx) = Xtez * w + c;
        yhat(testIdx) = score(testIdx) >= 0;
    end

    TP = sum((yhat == 1) & (y == 1));
    TN = sum((yhat == 0) & (y == 0));
    FP = sum((yhat == 1) & (y == 0));
    FN = sum((yhat == 0) & (y == 1));

    res = struct();
    res.tag = tag;
    res.N = numel(y);
    res.y = y;
    res.yhat = yhat;
    res.score = score;
    res.confusion = [TN FP; FN TP];
    res.accuracy = (TP + TN) / max(1, numel(y));
    res.sensitivity = TP / max(1, TP + FN);
    res.specificity = TN / max(1, TN + FP);
    res.precision = TP / max(1, TP + FP);
    res.f1 = 2 * res.precision * res.sensitivity / max(eps, res.precision + res.sensitivity);
    res.balancedAccuracy = 0.5 * (res.sensitivity + res.specificity);

    fprintf('\n%s\n', tag);
    fprintf('  N = %d\n', res.N);
    fprintf('  Accuracy           = %.3f\n', res.accuracy);
    fprintf('  Sensitivity        = %.3f\n', res.sensitivity);
    fprintf('  Specificity        = %.3f\n', res.specificity);
    fprintf('  Precision          = %.3f\n', res.precision);
    fprintf('  F1-score           = %.3f\n', res.f1);
    fprintf('  Balanced accuracy  = %.3f\n', res.balancedAccuracy);
end

function [w, c] = local_train_lda(X, y)
    X0 = X(y==0,:);
    X1 = X(y==1,:);

    mu0 = mean(X0,1).';
    mu1 = mean(X1,1).';

    S0 = cov(X0);
    S1 = cov(X1);

    if any(isnan(S0(:))) || rank(S0) < size(S0,1)
        S0 = S0 + 1e-6*eye(size(S0));
    end
    if any(isnan(S1(:))) || rank(S1) < size(S1,1)
        S1 = S1 + 1e-6*eye(size(S1));
    end

    n0 = size(X0,1);
    n1 = size(X1,1);
    Spooled = ((n0-1)*S0 + (n1-1)*S1) / max(1, n0+n1-2);
    Spooled = Spooled + 1e-6*eye(size(Spooled));

    invS = pinv(Spooled);
    w = invS * (mu1 - mu0);

    p0 = n0 / (n0+n1);
    p1 = n1 / (n0+n1);

    c = -0.5 * (mu1.' * invS * mu1 - mu0.' * invS * mu0) + log((p1+eps)/(p0+eps));
end
