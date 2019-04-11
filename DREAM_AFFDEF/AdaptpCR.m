function [pCR,lCR] = AdaptpCR(MCMCPar,CR,delta_tot,lCRold); 
% Updates the probabilities of the various crossover values

% Make CR to be a single vector
CR = CR(:);

% Determine lCR
for zz = 1:MCMCPar.nCR,
    % Determine how many times a particular CR value is used
    idx = find(CR==zz/MCMCPar.nCR);
    % This is used to weight delta_tot
    lCR(1,zz) = lCRold(1,zz) + size(idx,1);    
end;

% Adapt pCR using information from averaged normalized jumping distance
pCR = MCMCPar.seq * (delta_tot./lCR) / sum(delta_tot);

% Normalize pCR
pCR = pCR./sum(pCR);