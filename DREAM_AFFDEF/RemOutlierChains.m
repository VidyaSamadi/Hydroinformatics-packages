function [X,Sequences,hist_logp,outlier] = RemOutLierChains(X,Sequences,hist_logp,Iter,outlier,MCMCPar);
% Finds outlier chains and removes them when needed

% Determine the number of elements of L_density
[idx_end] = size(hist_logp,1); idx_start = floor(0.5*idx_end);

% Then determine the mean log density of the active chains
mean_hist_logp = mean(hist_logp(idx_start:idx_end,1:MCMCPar.seq));

% Initialize chain_id and Nid
chain_id = []; Nid = 0;

% Check whether any of these active chains are outlier chains
if strcmp(MCMCPar.outlierTest,'IQR_test'),
    % Derive the upper and lower quantile of the data
    Q1 = prctile(mean_hist_logp,75); Q3 = prctile(mean_hist_logp,25);
    % Derive the Inter quartile range
    IQR = Q1 - Q3;
    % Compute the upper range -- to detect outliers
    UpperRange = Q3 - 2 * IQR;
    % See whether there are any outlier chains
    chain_id = find(mean_hist_logp < UpperRange); Nid = size(chain_id,2);
end;

if strcmp(MCMCPar.outlierTest,'Grubbs_test'),
    % Test whether minimum log_density is outlier
    G = (mean(mean_hist_logp) - min(mean_hist_logp)) / std(mean_hist_logp);
    % Determine t-value of one-sided interval
    t2 = tinv(1 - 0.01/MCMCPar.seq,MCMCPar.seq-2)^2; % 95% interval
    % Determine the critical value
    Gcrit = ((MCMCPar.seq - 1)/sqrt(MCMCPar.seq)) * sqrt(t2/(MCMCPar.seq-2 + t2));
    % Then check this
    if (G > Gcrit), % Reject null-hypothesis
        % Indeed, an outlier chain
        chain_id = find(mean_hist_logp == min(mean_hist_logp)); Nid = 1;
    end;
end;

if strcmp(MCMCPar.outlierTest,'Mahal_test'),
    % Use the Mahalanobis distance to find outlier chains
    alpha = 0.01; UpperRange = ACR(MCMCPar.n,MCMCPar.seq-1,alpha);
    % Find which chain has minimum log_density
    [idx] = find(mean_hist_logp==min(mean_hist_logp)); idx = idx(1);
    % Set the other chains to 1
    ii = [1:MCMCPar.seq]; ii(idx) = 0; ii = find(ii>0);
    % Then check the Mahalanobis distance
    d1 = mahal(X(idx,1:MCMCPar.n),X(ii,1:MCMCPar.n));
    % Then see whether idx is an outlier in X
    if d1 > UpperRange,
        chain_id = idx; Nid = 1;
    end;
end;

if (Nid > 0),
    % Loop over each outlier chain
    for qq = 1:Nid,
        % Draw random other chain -- cannot be the same as current chain
        r_idx = find(mean_hist_logp==max(mean_hist_logp)); r_idx = r_idx(1);
        % Added -- update hist_logp -- chain will not be considered as an outlier chain then
        hist_logp(1:end,chain_id(qq)) = hist_logp(1:end,r_idx);
        % Jump outlier chain to r_idx -- Sequences
        Sequences(1,1:MCMCPar.n+2,chain_id(qq)) = X(r_idx,1:MCMCPar.n+2);
        % Jump outlier chain to r_idx -- X
        X(chain_id(qq),1:MCMCPar.n+2) = X(r_idx,1:MCMCPar.n+2);
        % Add to chainoutlier
        outlier = [outlier ; Iter chain_id(qq)]
    end;
else
    % Do nothing
end