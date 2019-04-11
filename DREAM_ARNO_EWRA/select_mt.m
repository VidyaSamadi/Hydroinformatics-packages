function [idx_sel] = select_mt(MCMCPar,log_p_y,log_pr1,log_sn);
% Pick one of the multi-try proposals in each chain based on the individual weights

% Work in log space
log_w = reshape(log_pr1,MCMCPar.mt,MCMCPar.seq) + reshape(log_p_y,MCMCPar.mt,MCMCPar.seq) + reshape(log_sn,MCMCPar.mt,MCMCPar.seq);

% Now for each column remove largest number
w = exp ( log_w - repmat(max(log_w),MCMCPar.mt,1) );

% Now select the best among the MCMCPar.mt proposals in each chain
for zz = 1:MCMCPar.seq,

    % Normalize each different column of w
    p = w(:,zz)./sum(w(:,zz)); 
    
    % Now select the "best" proposal from the normalized weights
    idx = randsample ( [ 1 : MCMCPar.mt ] , 1 , 'true' , p );
    
    % Now determine which element of y_mt this is for each chain
    idx_sel(zz,1) = idx + (zz-1) * MCMCPar.mt; 

end;