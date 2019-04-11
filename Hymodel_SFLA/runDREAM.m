% ----------------- DiffeRential Evolution Adaptive Metropolis algorithm -----------------------%
%                                                                                               %
% DREAM runs multiple different chains simultaneously for global exploration, and automatically %
% tunes the scale and orientation of the proposal distribution using differential evolution.    %
% The algorithm maintains detailed balance and ergodicity and works well and efficient for a    %
% large range of problems, especially in the presence of high-dimensionality and                %
% multimodality. 							                                                    %                     
%                                                                                               %
% DREAM developed by Jasper A. Vrugt and Cajo ter Braak                                         %
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%      input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain     %
%      Monte Carlo simulation, Water Resources Research, 44, W00B09, doi:10.1029/2007WR006720,  %
%      2008.                                                                                    %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by differential evolution with         %
%       self-adaptive randomized subspace sampling, International Journal of Nonlinear          %
%       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                                %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson, Equifinality of formal        %
%       (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling?, Stochastic     %
%       Environmental Research and Risk Assessment, 1-16, doi:10.1007/s00477-008-0274-y, 2009,  %
%       In Press.                                                                               %   
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution          %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
%   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater %
%       and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008.            %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, and J.M. Hyman, Differential evolution adaptive Metropolis   %
%       with snooker update and sampling from past states, SIAM journal on Optimization, 2009.  %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, and J.M. Hyman, Parallel Markov chain Monte Carlo simulation %
%       on distributed computing networks using multi-try Metropolis with sampling from past    %
%       states, SIAM journal on Scientific Computing, 2009.                                     %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
% All rights reserved.                                                                          %
%                                                                                               %
% Copyright 2008. Los Alamos National Security, LLC. This software was produced under U.S.      %
% Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is     %
% operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S.     %
% Government has rights to use, reproduce, and distribute this software.                        %
%                                                                                               %
% NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES A NY WARRANTY, EXPRESS OR  %
% IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to   %
% produce derivative works, such modified software should be clearly marked, so as not to       %
% confuse it with the version available from LANL.                                              %
%                                                                                               %
% Additionally, redistribution and use in source and binary forms, with or without              %
% modification, are permitted provided that the following conditions are met:                   %
% • Redistributions of source code must retain the above copyright notice, this list of         %
%   conditions and the following disclaimer.                                                    %
% • Redistributions in binary form must reproduce the above copyright notice, this list of      %
%   conditions and the following disclaimer in the documentation and/or other materials         %
%   provided with the distribution.                                                             %
% • Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL %
%   the U.S. Government, nor the names of its contributors may be used to endorse or promote    %
%   products derived from this software without specific prior written permission.              %
%                                                                                               %
% THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND   %
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES      %
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS %
% ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, %
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF   %
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)        %
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT %
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,       %
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                            %
%                                                                                               %
% MATLAB code written by Jasper A. Vrugt, Center for NonLinear Studies (CNLS)                   %
%                                                                                               %
% Written by Jasper A. Vrugt: vrugt@lanl.gov                                                    %
%                                                                                               %
% Version 0.5: June 2008                                                                        %
% Version 1.0: October 2008         Adaption updated and generalized CR implementation          %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% ----------------------------------------------------------------------------------------------%
%                                                                                               %
%           Note: DE-MC of ter Braak et al. (2006) is a special variant of DREAM                %
%                     It can be executed using the following options                            %
%                                                                                               %
%       MCMCPar.seq = 2 * MCMCPar.n;                                                            %
%       MCMCPar.DEpairs = 1;                                                                    %
%       MCMCPar.nCR = 1;                                                                        %
%       MCMCPar.eps = 0;                                                                        %
%       MCMCPar.outlierTest = 'None';                                                           % 
%       Extra.pCR = 'Fixed'                                                                     %
%                                                                                               %
% ----------------------------------------------------------------------------------------------%

% ----------------------------------------------------------------------------------------------%
%                                                                                               %
%   After termination of DREAM you can generate a 2D matrix with samples using the command:     %
%                                                                                               %
%       ParSet = GenParSet(Sequences,MCMCPar); if Extra.save_in_memory = 'Yes'                  %
%       otherwise                                                                               %
%       ParSet = GenParSet(Reduced_Seq,MCMCPar); if Extra.reduced_sample_collection = 'Yes'     %
%                                                                                               %
%       output.R_stat   The Gelman Rubin convergence diagnostic                                 %
%       output.AR       The acceptance percentage of candidate points                           %
%       output.CR       Optimized probabilities for individual crossover rates                  %
%       output.outlier  Outlier chain numbers                                                   %
%       hist_logp       Log density of last 50% of samples in each chain                        %
%       X               Final position of chains (population)                                   %
%                                                                                               %
% ----------------------------------------------------------------------------------------------%

% Different test examples
% example 1: n-dimensional banana shaped Gaussian distribution
% example 2: n-dimensional Gaussian distribution
% example 3: n-dimensional multimodal mixture distribution
% example 4: real-world example using hymod model

% Which example to run?
example = 4;

if example == 1, % n-dimensional banana shaped Gaussian distribution

    MCMCPar.n = 10;                         % Dimension of the problem (Nr. parameters to be optimized in the model)
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 3;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.ndraw = 100000;                 % Maximum number of function evaluations
    MCMCPar.steps = 10;                     % Number of steps
    MCMCPar.eps = 5e-2;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'Yes';% Thinned sample collection?
    Extra.T = 10;                           % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % Define the specific properties of the banana function
    Extra.mu   = [zeros(1,MCMCPar.n)];                      % Center of the banana function
    Extra.cmat = eye(MCMCPar.n); Extra.cmat(1,1) = 100;     % Target covariance
    Extra.imat = inv(Extra.cmat);                           % Inverse of target covariance
    Extra.bpar = [0.1];                                     % "bananity" of the target, see bananafun.m

    % What type of initial sampling
    Extra.InitPopulation = 'COV_BASED';
    % Provide information to do alternative sampling
    Extra.muX = Extra.mu;                                   % Provide mean of initial sample
    Extra.qcov = eye(MCMCPar.n) * 5;                        % Initial covariance
    % Save all information in memory or not?
    Extra.save_in_memory = 'No';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-Inf * ones(1,MCMCPar.n)]; ParRange.maxn = [Inf * ones(1,MCMCPar.n)];

    % Define the boundary handling
    Extra.BoundHandling = 'None';
    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'Banshp';
    % Define likelihood function
    option = 4;

end;

if example == 2,    % n-dimensional Gaussian distribution

    MCMCPar.n = 100;                        % Dimension of the problem
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.ndraw = 1000000;                % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 3;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 5e-2;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'Yes';% Thinned sample collection?
    Extra.T = 10;                           % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------
    
    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];
    %ParRange.minn = [9.9 * ones(1,MCMCPar.n)]; ParRange.maxn = [10 * ones(1,MCMCPar.n)];

    % Define the boundary handling
    Extra.BoundHandling = 'None';

    % ---------------------- Define covariance matrix ---------------------
    % Construct the d x d covariance matrix
    A = 0.5*eye(MCMCPar.n) + 0.5*ones(MCMCPar.n);
    % Rescale to variance-covariance matrix of interest
    for i=1:MCMCPar.n
        for j=1:MCMCPar.n
            C(i,j) = A(i,j)*sqrt(i*j);
        end
    end
    % Set to Extra
    Extra.qcov = C; Extra.muX = zeros(1,MCMCPar.n); Extra.invC = inv(C);
    % ---------------------------------------------------------------------

    Extra.save_in_memory = 'No';
    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define center of target
    Extra.mu = zeros(1,MCMCPar.n);          
    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Define modelName
    ModelName = 'normalfunc';
    % Define likelihood function
    option = 4;

end;

if example == 3,    % n-dimensional multimodal mixture distribution

    MCMCPar.n = 10;                         % Dimension of the problem
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.ndraw = 1000000;                % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 3;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 5e-2;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'Yes';% Thinned sample collection?
    Extra.T = 10;                           % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'COV_BASED';
    % Provide information to do alternative sampling
    Extra.muX = zeros(1,MCMCPar.n);         % Provide mean of initial sample
    Extra.qcov = eye(MCMCPar.n);            % Initial covariance

    Extra.Lam = eye(MCMCPar.n);             % covariance
    Extra.mu1 = -5 * ones(1,MCMCPar.n);     % center point of first density
    Extra.mu2 =  5 * ones(1,MCMCPar.n);     % center point of second density
    Extra.sigma = eye(MCMCPar.n);           % define the std of the target

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-Inf * ones(1,MCMCPar.n)]; ParRange.maxn = [Inf * ones(1,MCMCPar.n)];
    % Define the boundary handling
    Extra.BoundHandling = 'None';
    % Save in memory or not
    Extra.save_in_memory = 'No';

    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'mixturemodel';
    % Define likelihood function
    option = 1;

end;

if example == 4,    % AFFDEF rainfall - runoff model

    MCMCPar.n = 9;                          %Dimension of  the problem
    MCMCPar.seq = 15;              % Number of Markov Chains / sequences
    MCMCPar.ndraw = 4;                  % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0.99;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 1;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 2e-1;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [0.3 300 0.05 0.05 0.5 0.001 20000 .05 .05]; ParRange.maxn = [1.5 100000 10 10 20 0.1 800000 0.9 0.7];
    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Load the Leaf River data
    load obs94567.txt;

    Extra.MaxT = size(obs94567,1); 

    % Define the PET, Measured Streamflow and Precipitation.
    %Extra.PET = bound(1:Extra.MaxT,5);
    %Extra.Precip = sum(bound(1:Extra.MaxT,6:9),2);

    % Define the measured streamflow data
    Measurement.MeasData = obs94567(1:Extra.MaxT,1); Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'AFFDEF';
    % Define likelihood function -- Sum of Squared Error
    option = 2; Measurement.Sigma = ones(Measurement.N,1);
end;

% Scale of Delayed Rejection if used in algorithm
Extra.DR = 'No'; Extra.DRscale = 10;

% Run the MCMC algorithm
[Sequences,Reduced_Seq,X,output,hist_logp,prediction] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option);
% For postprocessing of results --> check comments on top of file


% Now create the output bounds
ParSet = genparset(Sequences,MCMCPar);

% Find the best solution
idx = find(ParSet(:,end)==max(ParSet(:,end))); 

% Call model to generate simulated data
evalstr = ['ModPred = ',ModelName,'(ParSet(idx(1),1:MCMCPar.n),Extra);'];
eval(evalstr);

% Compute the RMSE
RMSE = sqrt ( sum ( ( ModPred - Measurement.MeasData).^2) / Measurement.N);

% Take the last 20% or so -- make sure these are posterior samples -- look
% at R-statistic in output.R_stat;
Pars = ParSet(0.8*size(ParSet(:,1)):end,1:end);

% How many samples?
Ntot = size(Pars,1); 

% Assign simulation variables -- much faster
Smod = zeros(Measurement.N,Ntot); Spar = Smod;

% Define evalstr
evalstr = ['Spar(:,j) = ',ModelName,'(Pars(j,1:MCMCPar.n),Extra);']; 

% Now create the parameter and total uncertainty
for j = 1:Ntot,
    j
    % Call model to generate simulated data
    eval(evalstr);
    % Add the model error
    Smod(:,j) = Spar(:,j) + normrnd(0,RMSE,Measurement.N,1);
end;

% Now sort to get 95% ranges
for j = 1:Measurement.N,
    a = sort(Spar(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    minpar95(j,1) = a(floor(0.025 * Ntot)); maxpar95(j,1) = a(floor(0.975 * Ntot));
    % Same with total uncertainty
    a = sort(Smod(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    mintot95(j,1) = a(floor(0.025 * Ntot)); maxtot95(j,1) = a(floor(0.975 * Ntot));
end;

% Now plot the results 
plot([1:Measurement.N],minpar95,[1:Measurement.N],maxpar95); % These are the upper and lower parameter uncertainty bounds
hold on; plot([1:Measurement.N],mintot95,'R',[1:Measurement.N],maxtot95,'G'); % These are the upper and lower total uncertainty bounds
% Now add the observations
plot([1:Measurement.N],Measurement.MeasData,'r+');    

% Now calculate percentage inside the 95% bound --> should be about 95%
Contained = 1 - length(find(Measurement.MeasData<mintot95 | Measurement.MeasData > maxtot95))/Measurement.N;

% Note the ranges do not look very good --> because of homoscedasticity. If
% the error is assumed to increase with flow level results look much better