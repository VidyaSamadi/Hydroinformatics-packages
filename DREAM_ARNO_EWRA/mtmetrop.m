function [accept,ii,idx] = mtmetrop(MCMCPar,log_p_y,log_p_xref,log_pr1,log_pr2,log_sn,log_sn2,idx_sel);
% Multi-try Metropolis rule for acceptance or rejection of proposals

% And initialize accept with zeros
accept = zeros(MCMCPar.seq,1);

% So we determine, ratio =  p(y) * p(y --> X) / p(Xref) * p(Xref --> X)
% We go from X to y (using Z), and then from y to Xref
% p(y --> X) is equal to ||y - z||^(n-1) = log_sn
% p(Xref --> X) is equal to p(Xref --> y) * p(y --> X)

% Now determine total log-probability of mt proposals
log_w1 = log_p_y + log_pr1 + log_sn;

% Now determine total log-probability of xref, first p(y --> X), with y
% selected from mt proposals
log_sn_y = NaN(MCMCPar.mt * MCMCPar.seq,1);

% Now determine log_sn_y
for zz = 1:MCMCPar.seq,
    % Find the appropriate indices
    i_start = (zz-1) * MCMCPar.mt + 1; i_end = zz * MCMCPar.mt;
    % Determine log_sn_y from log_sn and idx_sel
    log_sn_y(i_start:i_end,1) = log_sn(idx_sel(zz));
end;

% Determine log_w2
log_w2 = log_p_xref + log_pr2 + log_sn_y + log_sn2;

% Calculate the weights for x and y
for zz = 1:MCMCPar.seq,

    % Find the appropriate indices
    i_start = (zz-1) * MCMCPar.mt + 1; i_end = zz * MCMCPar.mt;

    % Scale log density with maximum to avoid problems with zero
    m1 = max(log_w1(i_start:i_end,1)); m2 = max(log_w2(i_start:i_end,1)); m = max(m1,m2);

    % Determine the weights of y
    w1 = max ( exp ( log_w1(i_start:i_end,1) - m ), 1e-299 );

    % Determine the weights of x
    w2 = max ( exp ( log_w2(i_start:i_end,1) - m ), 1e-299 );

    % Calculate the Metropolis ratio
    alfa(zz,1) = sum(w1) / sum(w2);

end;

% Generate uniform random numbers -- to decide whether to accept or not
Z = rand(MCMCPar.seq,1);

% Find which alfa's are greater than Z
idx = find(alfa > Z);

% Find the original index of y
ii = idx_sel(idx);

% And indicate that these chains have been accepted
accept(idx,1) = 1;