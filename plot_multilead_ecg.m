function plot_multilead_ecg(X, fs, leadNames, figName, tLimitSec)
%PLOT_MULTILEAD_ECG Plot multi-lead ECG in stacked subplots.
%
% X can be leads x samples or samples x leads.

    if nargin < 5 || isempty(tLimitSec)
        tLimitSec = size(X,2) / fs;
    end

    if size(X,1) > size(X,2)
        X = X.';  % convert to leads x samples
    end

    nLeads = size(X,1);
    nShow = min(size(X,2), round(tLimitSec * fs));
    t = (0:nShow-1) / fs;

    figure('Name', figName, 'Color', 'w');
    for k = 1:nLeads
        subplot(nLeads,1,k);
        plot(t, X(k,1:nShow), 'k');
        grid on;
        if nargin >= 3 && ~isempty(leadNames)
            ylabel(leadNames{k});
        else
            ylabel(sprintf('L%d', k));
        end
        if k == 1
            title(figName);
        end
        if k == nLeads
            xlabel('Time (s)');
        else
            set(gca, 'XTickLabel', []);
        end
    end
end
