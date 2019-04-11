function [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);
% ------------- MT-DREAM with sampling from past and snooker updates: ------------------------- %
%                                                                                               %
% The code presented herein is a Markov Chain Monte Carlo algorithm that runs multiple chains   %
% in parallel for efficient posterior exploration. The algorithm, is based on the DREAM_(ZS)    %
% sampling scheme, but uses multi-try proposal sampling to increase sampling efficiency.        %
% Theory and numerical examples have been presented in Laloy and Vrugt (2012). The MT-DREAM_(ZS)%
% code is part of the DREAM MCMC package developed by Vrugt and coworkers.                      %
%                                                                                               %
% MT-DREAM_(ZS) developed by Jasper A. Vrugt and Eric Laloy                                     %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %
%                                                                                               %
% SYNOPSIS: [chain,X,Z,output] = mtdream_zs(MCMCPar,Func_name)                                  %
%           [chain,X,Z,output] = mtdream_zs(MCMCPar,Func_name,Extra)                            %
%           [chain,X,Z,output] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange)                   %
%           [chain,X,Z,output] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement)       %
%                                                                                               %
% Input:    MCMCPar = structure with DREAM parameters                                           %
%           Func_name = name of the function                                                    %
%           Extra = optional structure with arguments to be passed to function                  %
%           ParRange = optional structure with parameter ranges                                 %
%           Measurement = optional structure with measurement information                       %
%                                                                                               %
% Output:   chain = 3D array with Markov chain evolution                                        %
%           X = final position of chains and correponding density values                        %
%           Z = matrix with thinned sample history                                              %
%           output = structure with convergence properties, acceptance rate, CR values, etc.    %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   Laloy, E., and J.A. Vrugt, High-dimensional posterior exploration of hydrologic models      %
%       using multiple-try DREAM_(ZS) and high-performance computing, Water Resources Research, %
%       48, W01526, doi:10.1029/2011WR010608, 2012.                                             %
%                                                                                               %
%   ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater %
%       and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008             %
%                                                                                               %
%   Vrugt, J.A., E. Laloy, and C.J.F. ter Braak, DiffeRential Evolution Adaptive Metropolis     %
%       with Sampling from the Past and Subspace Updating, SIAM journal on Optimization         %
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Vrugt J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution           %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
%   ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%       input uncertainty in hydrologic modeling: Doing hydrology backward using Markov         %
%       chain Monte Carlo, Water Resour. Res., 44, W00B09, doi:10.1029/2007WR006720, 2008.      %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson, Equifinality of formal        %
%       (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling?, Stochastic     %
%       Environmental Research and Risk Assessment}, 23(7), 1011-1026,                          %
%       doi:10.1007/s00477-008-0274-y, 2009.                                                    %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by differential evolution with         %
%       self-adaptive randomized subspace sampling, International Journal of Nonlinear Sciences %
%       and Numerical Simulation}, 10(3), 273-290, 2009.                                        %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%     Copyright (C) 2011-2012  the authors                                                      %
%                                                                                               %
%     This program is free software: you can modify it under the terms of the GNU General       %
%     Public License as published by the Free Software Foundation, either version 3 of the      %
%     License, or (at your option) any later version.                                           %
%                                                                                               %
%     This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; %
%     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. %
%     See the GNU General Public License for more details.                                      %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
% Written by Jasper A. Vrugt: jasper@uci.edu                                                    %
%                                                                                               %
% Version 0.5: January 2009                                                                     %
% Version 1.0: April 2012         Maintenance update, postprocessing included                   %
% Version 1.1: January 2013       Approximate Bayesian Computation and simulation storage       %
% Version 1.2: September 2013     Added inference of measurement error as new option            %
% Version 1.3: May 2014           Parallellization using parfor (done if CPU > 1)               %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% --------------------------------------------------------------------------------------------- %
%                                                                                               %
%   Please check the directory /postprocessing which contains various                           %
%   script to postprocess the results                                                           %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% Calculate MCMCPar.steps and MCMCPar.m0
MCMCPar.steps = floor(MCMCPar.steps(MCMCPar)); MCMCPar.m0 = MCMCPar.m0(MCMCPar);

% Check how many input variables
if nargin < 5,
    % Define Measurement
    Measurement.MeasData = [];
end;

if nargin < 4,
    % Specify very large initial parameter ranges (minimum and maximum values)
    ParRange.minn = [ -Inf * ones(1,MCMCPar.n) ]; ParRange.maxn = [ Inf * ones(1,MCMCPar.n) ];
end;

if nargin < 3,
    % Define structure Extra to be empty
    Extra = [];
end;

% Now check how the measurement sigma is arranged (estimated or defined)
Measurement = CheckSigma(Measurement,MCMCPar);

% Check for simple setup errors
[stop,fid] = mtdream_zs_check(MCMCPar,Measurement);

% Check whether to stop the code due to setup errors?
if strcmp(lower(stop),'yes'),
    % Set output variables to be empty
    chain = []; X = []; output = []; Z = []; fx = [];
    % And return to main program
    return;
end;

% Determine updates
sampling_method = {'Parallel_Direction_Update','Snooker_Update'};

% Setup sequential or parallel environment
[MCMCPar] = mtdream_zs_calc_setup(MCMCPar);

% Initializes the main variables used in DREAM_ZS
[MCMCPar,chain,X,output,fx,fid_Z,delta_tot,pCR,lCR,CR,Iter,iteration,...
    T,T2,Table_JumpRate,iloc] = Initialize_mtdream_zs(MCMCPar,Func_name,Extra,...
    ParRange,Measurement);

% Initialize waitbar
h = waitbar(0,'Running MT-DREAM_{(ZS)} - Please wait...');

% Start timer
tic;

% Move prior population to posterior population ...
while (Iter < MCMCPar.ndraw),
    
    % Initialize totaccept;
    totaccept = 0;
    
    % Loop a number of times before calculating convergence diagnostic, etc.
    for gen_number = 1:MCMCPar.steps,
        
        % Initialize iteration
        T = T + 1; T2 = T2 + 1;
        
        % Define the current locations and associated posterior densities
        xold = X( 1 : MCMCPar.seq , 1 : MCMCPar.n); p_xold = X(1 : MCMCPar.seq , MCMCPar.n + 1); log_p_xold = ...
            X(1 : MCMCPar.seq , MCMCPar.n + 2);
        
        % Determine to do parallel direction or snooker update
        Update = sampling_method ( find ( mnrnd ( 1 , [ MCMCPar.parallelUpdate  1 - MCMCPar.parallelUpdate ] ) == 1 ) );
        
        % Genarate k multi-try proposals in each individual chain
        for uu = 1 : MCMCPar.mt,
            
            % Without replacement draw rows from Z for proposal creation
            R = randsample(MCMCPar.m, 2 * MCMCPar.DEpairs * MCMCPar.seq); Z = [];
            
            % Now read from file
            for zz = 1 : 2 * MCMCPar.DEpairs * MCMCPar.seq,
                % Load rows from "Z.bin"
                status = fseek(fid_Z,(R(zz)-1) * (MCMCPar.n + 2) * 8,'bof'); Z(zz,1:MCMCPar.n) = fread(fid_Z,MCMCPar.n,'double')';
            end;
            
            % Define idx
            idx = [uu : MCMCPar.mt : MCMCPar.mt * MCMCPar.seq];
            
            % Create the multi-try proposals
            [y_mt(idx,1:MCMCPar.n),CR_par(idx,1),log_sn(idx,1)] = offde(xold,Z,CR(:,gen_number),MCMCPar,Update,...
                Table_JumpRate,ParRange);
            
        end;
        
        % If an informative prior is used, prior density of each proposal of y_mt needs to be calculated
        [log_pr] = prior_pdf(MCMCPar,y_mt( 1 : MCMCPar.mt * MCMCPar.seq, 1:MCMCPar.n ) );
        
        % Now evaluate the model ( = pdf ) and return fx
        [fx_new] = Evaluate_model(y_mt,MCMCPar,Measurement,Func_name,Extra);
        
        % Calculate the likelihood and log-likelihood from fx
        [p_y,log_p_y] = CalcLikelihood(y_mt,fx_new,MCMCPar,Measurement);
        
        % Select one proposal in each chain using probability proportional to weight and set idx_tot
        [idx_sel] = select_mt(MCMCPar,log_p_y,log_pr,log_sn); idx_tot = zeros( MCMCPar.mt * MCMCPar.seq , 1 );
        
        % Modify CR -- to make sure that proposal and reference set use similar crossover values
        CR(:,gen_number) = CR_par(idx_sel,1); % --> theoretically needed, but statistically ok
        
        % Produce a reference set, x_ref, centered around the selected y
        for uu = 1 : MCMCPar.mt - 1,
            
            % Without replacement draw rows from Z for proposal creation
            R = randsample(MCMCPar.m, 2 * MCMCPar.DEpairs * MCMCPar.seq); Z = []; %Zoff = Z(R,1:MCMCPar.n);
            
            % Now read from file
            for zz = 1 : 2 * MCMCPar.DEpairs * MCMCPar.seq,
                % Load rows from "Z.bin"
                status = fseek(fid_Z,(R(zz)-1) * (MCMCPar.n + 2) * 8,'bof'); Z(zz,1:MCMCPar.n) = fread(fid_Z,MCMCPar.n,'double')';
            end;
            
            % Create idx and store which points have been created
            idx = [uu : MCMCPar.mt : MCMCPar.mt * MCMCPar.seq]; idx_tot(idx) = 1;
            
            % Create the reference set
            [xref(idx,1:MCMCPar.n),dummy,log_sn2(idx,1)] = offde(y_mt(idx_sel,1:MCMCPar.n),Z,CR(:,gen_number),...
                MCMCPar,Update,Table_JumpRate,ParRange);
            
        end;
        
        % Evaluate only the MCMCPar.mt - 1 proposals in each chain
        idx = find(idx_tot == 1); [fx_ref] = Evaluate_model(xref(idx,1:MCMCPar.n),MCMCPar,Measurement,Func_name,Extra);
        
        % Calculate the likelihood and log-likelihood from fx
        [p_xref(idx,1),log_p_xref(idx,1)] = CalcLikelihood(xref(idx,1:MCMCPar.n),fx_ref,MCMCPar,Measurement);
        
        % Augment xref, p_xref, and log_p_xref with current position of individual chains
        idx = find(idx_tot == 0); xref(idx,1:MCMCPar.n) = xold; p_xref(idx,1) = p_xold; log_p_xref(idx,1) = log_p_xold;
        log_sn2(idx,1) = 0;
        
        % If an informative prior is used the prior density of xref needs to be evaluated
        [log_pr2] = prior_pdf(MCMCPar,xref);
        
        % Apply the MTM-Metropolis rule (extended by JAV for generalized implementations with priors)
        [accept,ii_X,idx_X] = mtmetrop(MCMCPar,log_p_y,log_p_xref,log_pr,log_pr2,log_sn,log_sn2,idx_sel);
        
        % And update X and the model simulation
        X(idx_X,1:MCMCPar.n+2) = [y_mt(ii_X,1:MCMCPar.n) p_y(ii_X) log_p_y(ii_X)]; fx(:,idx_X) = fx_new(:,ii_X);
        
        % Check whether to add the current points to the chains or not?
        if (T == MCMCPar.T),
            
            % Store the current sample in chain
            iloc = iloc + 1; chain(iloc,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(X',1,MCMCPar.n+2,MCMCPar.seq);
            
            % Store the model simulations (if appropriate)
            [ T ] = mtdream_zs_store_results ( MCMCPar , fx , Measurement , 'a+' );
            
        end
        
        % Compute squared jumping distance for each CR value
        if strcmp(lower(MCMCPar.pCR),'yes');
            % Calculate the standard deviation of each dimension of X
            r = repmat(std(X(1:MCMCPar.seq,1:MCMCPar.n)),MCMCPar.seq,1);
            % Compute the Euclidean distance between new X and old X
            delta_normX = sum(((xold(1:end,1:MCMCPar.n) - X(1:end,1:MCMCPar.n))./r).^2,2);
            % Use this information to update sum_p2 to update N_CR
            [delta_tot] = CalcDelta(MCMCPar,delta_tot,delta_normX,CR(1:MCMCPar.seq,gen_number));
        end;
        
        % Check whether to append X to Z
        if (T2 == MCMCPar.k),
            % Append X to Z
            fclose(fid_Z); fid_Z = fopen('Z.bin','a+','n'); fwrite(fid_Z,X','double'); fclose(fid_Z); fid_Z = fopen('Z.bin','r');
            % Update MCMCPar.m
            MCMCPar.m = MCMCPar.m + MCMCPar.seq;
            % And set the T to 0
            T2 = 0;
        end;
        
        % How many candidate points have been accepted -- for Acceptance Rate
        totaccept = totaccept + sum(accept);
        
        % Update iteration
        Iter = Iter + MCMCPar.seq * ( 2 * MCMCPar.mt - 1);
        
        % Update the waitbar
        waitbar(Iter/MCMCPar.ndraw,h);
        
    end;
    
    % Reduce MCMCPar.steps to get rounded iteration numbers
    if (iteration == 2), MCMCPar.steps = MCMCPar.steps + 1; end;
    
    % Store Important Diagnostic information -- Acceptance Rate
    output.AR(iteration,1:2) = [ Iter 100 * totaccept/(MCMCPar.steps * MCMCPar.seq) ];
    
    % Store Important Diagnostic information -- Probability of individual crossover values
    output.CR(iteration,1:MCMCPar.nCR+1) = [ Iter pCR ];
    
    % Check whether to update individual pCR values
    if (Iter <= 0.1 * MCMCPar.ndraw),
        if strcmp(lower(MCMCPar.pCR),'yes'),
            % Update pCR values
            [pCR] = AdaptpCR(MCMCPar,delta_tot,lCR);
        end;
    end;
    
    % Generate CR values based on current pCR values
    [CR,lCRnew] = GenCR(MCMCPar,pCR); lCR = lCR + lCRnew;
    
    % Compute the R-statistic from chain using 50% burn-in
    [output.R_stat(iteration,1:MCMCPar.n+1)] = [ Iter Gelman(chain(max(1,floor(0.5*iloc)):iloc,1:MCMCPar.n,...
        1:MCMCPar.seq),MCMCPar) ];
    
    % Update the iteration
    iteration = iteration + 1;
    
    % Check whether we save or not?
    if strcmp(lower(MCMCPar.save),'yes');
        
        % Store in memory
        save MT-DREAM_ZS.mat
        
    end;
    
end;

% -------------------------------------------------------------------------

% Store CPU time
output.RunTime = toc;

% Variables have been pre-allocated --> need to remove zeros at end
[chain,Z,output,fx] = mtdream_zs_end(MCMCPar,Measurement,chain,output,iteration,iloc,pCR,fid_Z,nargout);

% Close the waitbar
close(h);
% Write final line of warning file
fprintf(fid,'--------- End of MT-DREAM_{(ZS)} warning file ---------\n');
% Close the warning file
fclose('all');
% Open the warning file
edit warning_file.txt
% Open the MCMC diagnostic file
edit MT-DREAM_ZS_diagnostics.txt

% -------------------------------------------------------------------------