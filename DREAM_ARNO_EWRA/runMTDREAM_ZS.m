% ----------------------------------------------------------------------------------------------%
%                                                                                               %
% DDDDDDDDDDDDDDD    RRRRRRRRRRRRRR     EEEEEEEEEEEEEEEE       AAAAA       MMM             MMM  %
% DDDDDDDDDDDDDDDD   RRRRRRRRRRRRRRR    EEEEEEEEEEEEEEEE       AAAAA       MMMM           MMMM  %
% DDD          DDD   RRR          RRR   EEE                   AAA AAA      MMMMM         MMMMM  %
% DDD          DDD   RRR          RRR   EEE                   AAA AAA      MMMMMM       MMMMMM  %
% DDD          DDD   RRR          RRR   EEE                  AAA   AAA     MMM MMM     MMM MMM  %
% DDD          DDD   RRR          RRR   EEE                  AAA   AAA     MMM  MMM   MMM  MMM  %
% DDD          DDD   RRRRRRRRRRRRRRRR   EEEEEEEEEEEEEEEE    AAA     AAA    MMM   MMM MMM   MMM  %
% DDD          DDD   RRRRRRRRRRRRRRRR   EEEEEEEEEEEEEEEE    AAAAAAAAAAA    MMM    MMMM     MMM  %
% DDD          DDD   RRR          RRR   EEE                AAA       AAA   MMM             MMM  %
% DDD          DDD   RRR          RRR   EEE                AAA       AAA   MMM             MMM  %
% DDD          DDD   RRR          RRR   EEE               AAA         AAA  MMM             MMM  %
% DDD          DDD   RRR          RRR   EEE               AAA         AAA  MMM             MMM  %
% DDDDDDDDDDDDDDDD   RRR          RRR   EEEEEEEEEEEEEEEE  AAA         AAA  MMM             MMM  %
% DDDDDDDDDDDDDDD    RRR          RRR   EEEEEEEEEEEEEEEE  AAA         AAA  MMM             MMM  %
%                                                                                               %
% ----------------------------------------------------------------------------------------------%

% ------------- MT-DREAM with sampling from past and snooker updates -------------------------- %
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
% The directory \PostProcessing contains a script "PostProcMCMC" that will compute various      %
% posterior statistics (MEAN, STD, MAP, CORR) and create create various plots including,        %
% marginal posterior parameter distributions, R_stat convergence diagnostic, two-dimensional    %
% correlation plots of the posterior parameter samples, chain convergence plots, and parameter  %
% and total posterior simulation uncertainty ranges (interval can be specified)                 %
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
% Version 1.0: April 2011         Maintenance update, explicit treatment of prior distribution  %
% Version 1.1: January 2013       Approximate Bayesian Computation and simulation storage       %
% Version 1.2: September 2013     Added inference of measurement error as new option            %
% Version 1.3: May 2014           Parallellization using parfor (done if CPU > 1)               %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% Different test examples
% example 1:    n-dimensional banana shaped Gaussian distribution
% example 2:    n-dimensional Gaussian distribution
% example 3:    n-dimensional multimodal mixture distribution
% example 4:    real-world example using hymod rainfall - runoff model (HYMOD code in MATLAB)
% example 5:    real-world example using hymod rainfall - runoff model (HYMOD code in FORTRAN)
% example 6:    rainfall-runoff model with generalized log-likelihood function
% example 7:    HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters
% example 8:    multivariate student t distribution
% example 9:    rainfall-runoff model with Whittle's likelihood function
% example 10:   real-world example using SAC-SMA rainfall - runoff model
% example 11:   one-dimensional toy example Approximate Bayesian Computation
% example 12:   ABC inference for hydrologic model
% example 13:   Three-dimensional soil moisture inference from GPR travel time data
% example 14:   Inference of crosshole ground penetrating radar slowness distribution based on a straight-ray approximation using the discrete cosine transform.
% example 15:   Hymod rainfall - runoff model with estimation of measurement error
    
% Most of the examples herein are used to illustrate how the code works.
% Yet, for some of these problems it is much better to run DREAM and DREAM_ZS as
% this will provide many more posterior samples and hence better posterior
% estimates. The multi-try sampling per definition stores far fewer samples
% for a given number of function evaluations, and hence this affects the
% posterior estimates. Hence please use DREAM and DREAM_ZS for
% simpler problems. For distributed problems storage is sufficient as the
% the number of function evaluations will typically be very large.

% Clear memory
clear all; clc;

% Which example to run?
example = 4;

global DREAM_dir EXAMPLE_dir CONV_dir

% Store working directory and subdirectory containing the files needed to run this example
DREAM_dir = pwd; EXAMPLE_dir = [pwd '\example_' num2str(example)]; CONV_dir = [pwd '\diagnostics'];

% Add subdirectory to search path
addpath(EXAMPLE_dir); addpath(CONV_dir); addpath(DREAM_dir);

% Recommended parameter settings
MCMCPar.seq = 3;                        % Number of Markov chains / chain (for high dimensional and highly nonlinear problems, larger values work beter!!)
MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
MCMCPar.nCR = 3;                        % Number of crossover values used
MCMCPar.k = 10;                         % Thinning parameter for appending X to Z (make sure that enough posterior samples are being stored!!)
MCMCPar.mt = 5;                         % Number of multi-try trials (try 3 as well!)
MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
MCMCPar.steps = inline('MCMCPar.ndraw/(250 * MCMCPar.seq)'); % Number of steps before calculating convergence diagnostics
MCMCPar.m0 = inline('10 * MCMCPar.n');  % Initial size of matrix Z
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes
MCMCPar.pCR = 'Yes';                    % Adaptive tuning of crossover values (Yes or No)
MCMCPar.Restart = 'No';                 % Restart run (Yes or No)
MCMCPar.modout = 'Yes';                 % Return model (function) simulations of samples Yes or No)?
MCMCPar.save = 'Yes';                   % Save output during the run (Yes or No)
MCMCPar.parallel = 'No';                % Parallel evaluation or not (Yes or No)
MCMCPar.ABC = 'No';                     % Approximate Bayesian Computation or Not?
MCMCPar.IO = 'No';                      % Input-output writing of model (Yes or No)

% NOTE -- MT-DREAM(ZS) is especially designed for parallel MCMC simulation.
% I recommend using DREAM and DREAM(ZS) for simple problems.

% -------------------------------------------------------------------------
%     Likelihood functions can be found in the function CompDensity.m
% Some of these functions (2 and 7) require the measurement error sigma to 
% be defined. One can set "Measurement.Sigma" to be a single value. This 
% assumes homoscedastic errors. Or one can assume heteroscedasticity by
% using a vector of similar length as the Measurement.MeasData. As an 
% alternative one can also define a model for the measurement error, for
% instance, Measurement.Sigma = inline('a * x + b'), where "a" and "b" are
% additional parameters to be estimated, and "x" is Measurement.MeasData.
% Example 15 gives an illustration how to do this in practice.
% -------------------------------------------------------------------------

if example == 1, % n-dimensional banana shaped Gaussian distribution

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,   %
    %       Accelerating Markov chain Monte Carlo simulation by differential evolution with     %
    %       self-adaptive randomized subspace sampling, International Journal of Nonlinear      %
    %       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                            %
    %                                                                                           %
    %   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution      %
    %       Metropolis algorithm for optimization and uncertainty assessment of hydrologic      %
    %       model parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003. %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Application specific parameter settings
    MCMCPar.n = 10;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'COV';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % Define the likelihood function --> log-density from model

    % (GIVEN THE NONLINEARIRY OF THE BANANA SHAPED CURVE --> THINK ABOUT INCREASING MCMCPAR.SEQ)

    % Define Func_name
    Func_name = 'Banana_func';

    % Provide information to do initial sampling ("COV")
    MCMCPar.mu = zeros(1,MCMCPar.n);        % Provide mean of initial sample
    MCMCPar.cov = 5 * eye(MCMCPar.n);       % Initial covariance

    % Define the specific properties of the banana function
    Extra.cov  = eye(MCMCPar.n); Extra.cov(1,1) = 100;     % Target covariance
    Extra.invC = inv(Extra.cov);                           % Inverse of target covariance

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra);

end;

if example == 2,    % n-dimensional Gaussian distribution

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,   %
    %       Accelerating Markov chain Monte Carlo simulation by differential evolution with     %
    %       self-adaptive randomized subspace sampling, International Journal of Nonlinear      %
    %       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                            %
    %                                                                                           %
    %   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker     %
    %       updater and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9,      %
    %       2008.                                                                               %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 100;                        % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 1000000;                % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % Model returns log-density

    % Define Func_name
    Func_name = 'normalfunc';

    % Give the parameter ranges (minimum and maximum values) (mostly overdispersed)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];

    % Give the parameter ranges (minimum and maximum values) (underdispersed)
    % ParRange.minn = [9.9 * ones(1,MCMCPar.n)]; ParRange.maxn = [10 * ones(1,MCMCPar.n)];

    % ------ Define covariance matrix of target distribution --------------

    % Construct the d x d covariance matrix
    A = 0.5 * eye(MCMCPar.n) + 0.5 * ones(MCMCPar.n);
    % Rescale to variance-covariance matrix of interest
    for i = 1:MCMCPar.n
        for j = 1:MCMCPar.n
            C(i,j) = A(i,j) * sqrt(i * j);
        end
    end

    % Store inverse of covariance for "normalfunc"
    Extra.invC = inv(C);

    % ---------------------------------------------------------------------

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange);

end;

if example == 3,    % n-dimensional multimodal mixture distribution

    % ---------------------------- Check the following two papers ----------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,   %
    %       Accelerating Markov chain Monte Carlo simulation by differential evolution with     %
    %       self-adaptive randomized subspace sampling, International Journal of Nonlinear      %
    %       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                            %
    %                                                                                           %
    %   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker     %
    %       updater and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9,      %
    %       2008.                                                                               %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 25;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 1500000;                % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'COV' ;                 % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 1;                        % Model returns directly the log-probability density
    MCMCPar.seq = 5;

    % Define Func_name
    Func_name = 'mixturemodel';

    % Provide information to do covariance sampling ("COV")
    MCMCPar.mu = zeros(1,MCMCPar.n);        % Provide mean of initial sample
    MCMCPar.cov = eye(MCMCPar.n);           % Initial covariance

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name);

end;

if example == 4,    % HYMOD rainfall - runoff model (coded in MATLAB)

    % ---------------------------- Check the following 3 papers ------------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of  %
    %      input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain %
    %      Monte Carlo simulation, Water Resources Research, 44, W00B09,                        %
    %      doi:10.1029/2007WR006720, 2008.                                                      %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson, Equifinality of formal    %
    %       (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling?, Stochastic %
    %       Environmental Research and Risk Assessment, 23(7), 1011-1026, 				        %
    %       doi:10.1007/s00477-008-0274-y, 2009                                                 %
    %                                                                                           %
    %   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution      %
    %       Metropolis algorithm for optimization and uncertainty assessment of hydrologic      %
    %       model parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003. %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 15;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 3;                        % Define likelihood function -- Sum of Squared Error

    % Define Func_name
    Func_name = 'ARNO';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1 10 10 1 100 1000 50 0 0 0.01 0 0 0 1.5 0 ]; ParRange.maxn = [3 1000 100 3 10000 100000 600 300 100 1 10 10 5 5 10];

    % Load the Leaf River data
    load Hisflow25.txt;
    load Hisrain25.txt;
    load Hispan25.txt;

    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 2557;

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.E = Hispan25(1:Extra.MaxT,1); Extra.P = Hisrain25(1:Extra.MaxT,1);

    % Area factor to translate HYMOD output in mm/d to m3/s (calibration data); (area Leaf River is 1944 km2)
    Extra.F = 2300 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);

    % Define the measured streamflow data
    Measurement.MeasData = Hisflow25(1:Extra.MaxT,1);

    % We need to specify the Measurement error of the data in Measurement.Sigma
    % With option 3, Measurement.Sigma is integrated out the likelihoon function
    % With any other option, Sigma needs to be defined

    % We can estimate the measurement error directly if we use temporal differencing
    % The function MeasError provides an estimate of error versus flow level
    % out = MeasError(Measurement.MeasData;
    % For the Leaf River watershed this results in a heteroscedastic error
    % that is about 10% of the actual measured discharge value, thus
    % You can check this by plotting out(:,1) versus out(:,2)
    % Measurement.Sigma = 0.1*Measurement.MeasData; % --> option 2

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);

end;

if example == 5,    % HYMOD rainfall - runoff model (coded in FORTRAN)

    % Problem specific parameter settings
    MCMCPar.n = 5;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 3;                        % Define likelihood function -- Sum of Squared Error
    MCMCPar.IO = 'Yes';                     % Input-output writing of model (Yes or No)
    
    % Define Func_name
    Func_name = 'hymodFORTRAN';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10]; ParRange.maxn = [500 2.00 0.99 0.10 0.99];

    % Leaf River data -- forcing conditions not needed --> externally loaded by FORTRAN executable
    load bound.txt;

    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795;

    % Define the measured streamflow data -- use to compute likelihood function
    Measurement.MeasData = bound(65:Extra.MaxT,4);

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);

end;

if example == 6,    % Rainfall-runoff model with generalized log-likelihood

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % G. Schoups, J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption of          %
    %     accuracy and efficiency of Markov Chain Monte Carlo simulation by inaccurate          %
    %     numerical implementation of conceptual hydrologic models, Water Resources             %
    %     Research, 46, W10530, doi:10.1029/2009WR008648.                                       %
    %                                                                                           %
    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and         %
    %     predictive inference of hydrologic models with correlated, heteroscedastic and        %
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933.  %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 11;                         % Dimension of the problem
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 8;                        % Generalized likelihood function

    % Define Func_name
    Func_name = 'hmodel';

    % Give the parameter ranges (minimum and maximum values)
    %parno:       1     2     3     4     5     6     7      8    9     10    11   12   13   14   15   16   17   18   19   20  21
    %parname:     fA    Imax  Smax  Qsmax alE   alS   alF    Pf   Kfast Kslow std0 std1 beta xi   mu1  phi1 phi2 phi3 phi4 K   lambda
    Extra.fpar = [1     0     100   10    100   1e-6  1e-6   0    2     70    0.1  0    0    1    0    0    0    0    0    0   1];
    % Minimum parameter values
    parmin =     [1     0     10    0     1e-6  1e-6 -10     0    0     0     0    0   -1    0.1  0    0    0    0    0    0   0.1 ];
    % maximum parameter values
    parmax =     [1     10    1000  100   100   1e-6  10     0    10    150   1    1    1    10   100  1    1    1    1    1   1];
    % Select the parameters to be sampled
    Extra.idx_vpar = [2 3 4 5 7 9 10 11 12 13 16];

    % And define the associated Parameter ranges
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.Precip = daily_data(:,4); Extra.Ep = daily_data(:,5);

    % Define the measured streamflow data
    Measurement.MeasData = daily_data(Extra.idx,6);

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);

end;

if example == 7,	% HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters

    % -------------------------------- Check the following paper ------------------------------ %
    %                                                                                           %
    %   B. Scharnagl, J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse          %
    %	modeling of soil water dynamics at the field scale: using prior information             %
    %	on soil hydraulic properties, Hydrology and Earth System Sciences.                      %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'PRIOR';                % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % The model directly returns the log-density
    MCMCPar.IO = 'Yes';                     % Input-output writing of model (Yes or No)
    
    % Define model name
    Func_name = 'HYDRUS';

    % Define feasible parameter space (minimum and maximum values)
    %				1		2		3				4			5			6		7
    %				[thetar	thetas	log10(alpha)	log10(n)	log10(Ks)	L		hLB
    ParRange.minn =	[0.0430	0.4090	-2.5528			0.1790		-2.2366		-5.4900	-250];
    ParRange.maxn =	[0.0910 0.4810	-2.0706			0.2670		-0.0800		6.2700	-50];

    % Provide observational data and data needed to modify the initial and boundary conditions
    [Extra] = ProvideData;

    % Specify the prior distributions for the various parameters
    MCMCPar.prior_marginal = {  'normrnd(0.0670,0.0060)',...
        'normrnd(0.4450,0.0090)',...
        'normrnd(-2.310,0.0600)',...
        'normrnd(0.2230,0.0110)',...
        'normrnd(-1.160,0.2700)',...
        'normrnd(0.3900,1.4700)',...
        'unifrnd(-250,-50)'};

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange);

end;    

if example == 8,    % multivariate student t distribution with 60 degrees of freedom

    % ---------------------------- Check the following paper ---------------------------------- %
    %                                                                                           %
    %   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker     %
    %       updater and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9,      %
    %       2008.                                                                               %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 25;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 10;                         % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 1;                        % Model returns density

    % Define Func_name
    Func_name = 'multi_student';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];

    % ---------------------- Define covariance matrix ---------------------

    % Construct the dxd correlation matrix
    Extra.corr = 0.5 * eye(MCMCPar.n) + 0.5 * ones(MCMCPar.n);

    % Define the example specific properties used to compute output
    Extra.mu = zeros(1,MCMCPar.n);

    % How many degrees of freedom of student distribution used as target function?
    Extra.df = 60;

    % Make sure C is a valid covariance matrix
    [Extra.R,err] = cholcov(Extra.corr,0);

    % Define logSqrtDetC
    Extra.logSqrtDetC = sum(log(diag(Extra.R)));

    % Define dimensionality
    Extra.d = MCMCPar.n;

    % ---------------------------------------------------------------------

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange);

end;

if example == 9,    % Rainfall-runoff model with Whittle's log-likelihood (spectral analysis)

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % G. Schoups, J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption of          %
    %     accuracy and efficiency of Markov Chain Monte Carlo simulation by inaccurate          %
    %     numerical implementation of conceptual hydrologic models, Water Resources             %
    %     Research, 46, W10530, doi:10.1029/2009WR008648.                                       %
    %                                                                                           %
    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and         %
    %     predictive inference of hydrologic models with correlated, heteroscedastic and        %
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933.  %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % and for Whittle's likelihood function and application:

    % ----------------------------------------------------------------------------------------- %
    %                                                                                           %
    % A. Montanari, and E. Toth (2007), Calibration of hydrological models in the spectral      %
    % domain: An opportunity for scarcely gauged basins?, Water Resources Research, 43, W05434, %
    % doi:10.1029/2006WR005184.                                                                 %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 9;                        % Whittle's likelood function -- spectral analysis of data

    % Define Func_name
    Func_name = 'hmodel';

    % Give the parameter ranges (minimum and maximum values)
    %parno:       1     2     3     4     5     6     7      8    9     10    11   12   13   14   15   16   17   18   19   20  21
    %parname:     fA    Imax  Smax  Qsmax alE   alS   alF    Pf   Kfast Kslow std0 std1 beta xi   mu1  phi1 phi2 phi3 phi4 K   lambda
    Extra.fpar = [1     0     100   10    100   1e-6  1e-6   0    2     70    0.1  0    0    1    0    0    0    0    0    0   1];
    parmin =     [1     0     10    0     1e-6  1e-6 -10     0    0     0     0    0   -1    0.1  0    0    0    0    0    0   0.1 ];
    parmax =     [1     10    1000  100   100   1e-6  10     0    10    150   1    1    1    10   100  1    1    1    1    1   1];
    Extra.idx_vpar = [2 3 4 5 7 9 10];
    % Define parameter ranges
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.Precip = daily_data(:,4); Extra.Ep = daily_data(:,5);

    % Define the measured streamflow data
    Measurement.MeasData = daily_data(Extra.idx,6);

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);

end;

if example == 10,    % SAC-SMA hydrological model

    % Problem specific parameter settings
    MCMCPar.n = 13;                         % Dimension of the problem
    MCMCPar.ndraw = 250000;                 % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 3;                        % log-density evaluation

    % Define Func_name
    Func_name = 'sacmodel';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 1.0 0.1 0.0 0.0 1.0 1 1.0 1.00 1.00 0.01 0.0001 0.0]; %0  ];
    ParRange.maxn = [150 150 0.5 0.1 0.4 250 5 500 1000 1000 0.25 0.0250 0.6]; %50 ];

    % Opt derived with AMALGAM algorithm
    x_opt = [16.667930	34.203931	0.308337	0.000046	0.257849	249.984817	3.045226	264.282574	27.084401	89.329993	0.249591	0.018338	0.099927	]; RMSE = 13.22;

    % Define data structures for use in computation of posterior density
    load bound.txt; Extra.Obs = bound(1:795,5:9);

    % Store in array Measurement.MeasData
    Measurement.MeasData = bound(65:795,4);

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);

end;

if example == 11,    % Simple 1D mixture distribution -- Approximate Bayesian Computation

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % M. Sadegh, and J.A. Vrugt (2013), Approximate Bayesian computation using DREAM: Theory,   %
    %     Numerical Implementation and Case Studies, Water Resources Research, In Prep.         %
    %                                                                                           %
    % B.M. Turner, and P.B. Sederberg (2013), Approximate Bayesian computation with             %
    %     differential evolution, Journal of Mathematical Psychology, In Press.                 %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 1;                          % Dimension of the problem
    MCMCPar.ndraw = 100000;                 % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 11;                       % ABC likelihood function (standardized)
    MCMCPar.ABC = 'Yes';                    % Specify that we perform ABC
    MCMCPar.delta = 0.025;                  % Delta of the noisy ABC implementation
    MCMCPar.rho = inline('X - Y');          % Define the distance - in this case the difference

    % Define Func_name
    Func_name = 'ABC_func';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-10]; ParRange.maxn = [10];

    % Define Measurement.MeasData --> "Y" in paper
    Measurement.MeasData = 0;

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,[],ParRange,Measurement);

end;

if example == 12, % Rainfall runoff modeling problem using Schoups and Vrugt model (2010)

    % ---------------------------- Check the following 3 papers ------------------------------- %
    %                                                                                           %
    % J.A. Vrugt, and M. Sadegh, Towards diagnostic model calibration and evaluation:           %
    %     Appproximate Bayesian computation, Water Resources Research, 2012, In Review.         %
    %                                                                                           %
    % M. Sadegh, and J.A. Vrugt (2013), Approximate Bayesian computation using DREAM: Theory,   %
    %     Numerical Implementation and Case Studies, Water Resources Research, In Prep.         %
    %                                                                                           %
    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and         %
    %     predictive inference of hydrologic models with correlated, heteroscedastic and        %
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933.  %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 10;                       % ABC likelihood function (standardized)
    MCMCPar.ABC = 'Yes';                    % Specify that we perform ABC
    MCMCPar.delta = 0.025;                  % Delta of the noisy ABC implementation
    MCMCPar.rho = inline('X - Y');          % Define the distance - in this case the difference

    % Give the parameter ranges (minimum and maximum values)
    %parno         1     2      3      4     5      6      7
    %parname:     Imax  Smax  Qsmax   alE   alF   Kfast  Kslow

    % Minimum parameter values
    parmin =     [ 0     10     0    1e-6   -10     0      0  ];
    % maximum parameter values
    parmax =     [ 10   1000   100   100     10     10    150 ];

    % Select the parameters to be sampled
    Extra.idx_vpar = [1:7];

    % And define the associated Parameter ranges
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.P = daily_data(:,4); Extra.Ep = daily_data(:,5);

    % Define the measured streamflow data
    MeasData = daily_data(Extra.idx,6);

    % Now calculate the summary metrics
    [Measurement.MeasData] = CalcMetrics(MeasData,Extra.P(Extra.idx))';

    % Define Func_name
    Func_name = 'rainfall_runoff';

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);

end;

if example == 13,   % Three-dimensional water tracer distribution in the vadose zone from GPR data

    % ---------------------------- Check the following paper ---------------------------------- %
    %                                                                                           %
    % E. Laloy, N. Linde, and J.A. Vrugt, Mass conservative three-dimensional water tracer      %
    %     distribution from Markov chain Monte Carlo inversion of time-lapse ground-penetrating %
    %     radar data, Water Resour. Res., 48, W07510, doi:10.1029/2011WR011238.                 %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.ndraw = 500000;                 % Maximum number of function evaluations
    MCMCPar.T = 10;                         % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % The model directly returns the log-density

    % Define Func_name
    Func_name = 'Legendre_GPR';

    % Problem dimensions
    Extra.order = 5 ;                       % Order of the legendre moments used for reconstructing the image
    Extra.nalpha = 15;                      % Number of alpha parameters considered (problem-specific),
    % for instance, Legendre moments of order 4 correspond herein to 8 alpha parameters,
    % order 5 to 15, order 6 to 24, and order 7 to 35 (see calculation in init_ConsSVD_Tau.m)
    Extra.nbeta = 4 ;                       % Four shape parameters in 2D
    Extra.error = 0.5;                      % Standard deviation of Gaussian traveltime error (ns)
    MCMCPar.n = Extra.nalpha + Extra.nbeta; % Dimension of the problem (number of parameters to be estimated)

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn=[zeros(1,Extra.nbeta) , zeros(1,Extra.nalpha) - 0.1];
    ParRange.maxn=[zeros(1,Extra.nbeta) + 14 , zeros(1,Extra.nalpha) + 0.1];

    % The x-axis is varying the fastest
    x=0.0:0.1:3;                            % Boundaries of uniform x-grid (3 by 3 m grid); 0.1 m discretization
    z=0:0.1:3;                              % Boundaries of uniform z-grid

    sourcex = 0.01;                         % x-position of source
    sourcez = 0.05 : 0.1 : 3;               % z-positions of sources
    receiverx = 2.99;                       % x-position of receivers
    receiverz = 0.05 : .1 : 3;              % z-positions of receivers
    nsource = length(sourcez); nreceiver = length(receiverz);

    % Calculate acquisition geometry (multiple-offset gather)
    for j=1:nsource
        for i=1:nreceiver
            data((j-1)*nreceiver+i,:)=[sourcex sourcez(j) receiverx receiverz(i)];
        end
    end

    % Calculate forward modeling kernel, Courtesy Dr. James Irving, UNIL
    Extra.J = tomokernel_straight(data,x,z); % Distance of ray-segment in each cell for each ray

    % Prepare moisture plume and grid parameters
    Extra.dim=2;                            %2D model
    dx = 0.1;                               % Discretization of original cells (i.e., to be used in forward modeling)
    dz = 0.1;
    nz = 30;                                % Original discretization of water saturation model
    nx = 40;
    wcon = load('Sw.dat');                  % Load water saturation model: Courtesy of Dr. Michael Kowalsky (LBNL), Kowalsky et al. (2005; WRR)
    Extra.por = 0.36;                       % Porosity field
    for j = 1 : nz
        for i = 6 : nx - 5                  % Make model 3 by 3 m
            wcont(j,i-5) = wcon(nz-j+1,i);
        end
    end
    nx = 30;                                % new nx
    wcont = wcont * Extra.por;              % Transform saturation into water content

    % Compute grid parameters
    Extra.GridPar = [nx,nz,dx,dz];
    Extra.xiv = [0 + dx/2 : dx : nx*dx-dx/2];
    Extra.ziv = [0 + dz/2 : dz : nz*dz-dz/2];
    [Extra.xi,Extra.zi] = meshgrid(Extra.xiv,Extra.ziv);

    % Create plume background and increment
    Extra.W0 = repmat(wcont(:,1),1,nx);
    Extra.dW = wcont - Extra.W0 ;
    % Make sure the increment is always >=0
    Extra.dW(Extra.dW<0) = 0;
    % And build the true plume
    Extra.W1 = Extra.W0 + Extra.dW;
    % Compute added mass (m3)
    Extra.mass = sum((Extra.dW(:))) * dx * dz;

    % Translate water content model into slowness using the Refractive Index/CRIM model
    Extra.pw = 81;                          % Permittivity of water
    Extra.pa = 1;                           % Permittivity of air
    Extra.ps = 5;                           % Permittivity of mineral grains

    % Create slowness model (ns/m)
    slowtrue = Extra.W1 * sqrt(Extra.pw) + (Extra.por-Extra.W1) * sqrt(Extra.pa) + (1-Extra.por) * sqrt(Extra.ps); % sqrt of effective permittivity
    slowtrue = slowtrue / 0.3;              % True slowness field

    % Simulate data for true slowness model
    % Add Gaussian uncorrelated noise with a standard deviation of 0.5 ns.
    Extra.datasim = Extra.J * slowtrue(:) + Extra.error * randn(nsource*nreceiver,1); % Simulate data

    % Prepare matrix of prior constraints, the vectors and matrices needed
    % for its SVD decomposition, and the Tau matrix to be used for plume reconstruction from its inferred Legendre moments
    [Extra] = init_ConsSVD_Tau(Extra) ;
    % ---------------------------------------------------------------------

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange);

end;

if example == 14,   % Crosshole GPR slowness distribution based on a straight-ray approximation using the discrete cosine transform. 
                    % The problem is simplified compared to the paper cited below as it considers a problem in which the true model 
                    % represent a model with the same dimension as the inverse parameterization and it uses straight rays.

    % ### Results can be visualized by the function visualize_results.m ###

    % ---------------------------- Check the following paper ---------------------------------- %
    %                                                                                           %
    %  N. Linde, and J.A. Vrugt, Spatially distributed soil moisture from traveltime            %
    %       observations of crosshole ground penetrating radar using Markov chain Monte Carlo,  %
    %       Vadose Zone Journal, 2013                                                           %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %
       
    Extra.parx = 8;                         % Inversion parameters in x-direction (DCT order)
    Extra.parz = 8;                         % Inversion parameters in z-direction (DCT order)
                                            % 64 dimensions in total
    
    % Problem specific parameter settings
    MCMCPar.n = Extra.parx * Extra.parz;    % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 300000;                 % Maximum number of function evaluations
    MCMCPar.T = 5;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'COV';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % The model directly returns the log-density
    %MCMCPar.save='Yes';

    % Define Func_name
    Func_name = 'DCT_GPR';
    Extra.error = 0.5;                      % Standard deviation of Gaussian traveltime error

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn(1) = 30 * 0.7696; ParRange.maxn(1) = 30 * 1.301;              % Corresponds to 1/0.05 to 1/0.17 ns/m in logarithmic units

    % Set-up forward kernel
    % The x-axis is varying the fastest
    x = 0 : 0.1 : 3;                        % Boundaries of uniform x-grid (3 by 3 m grid); 0.1 m discretization
    z = 0 : 0.1:3;                          % Boundaries of uniform z-grid

    sourcex = 0.01;                         % x-position of source
    sourcez = 0.05 : 0.1 : 3.;              % z-positions of sources
    receiverx = 2.99;                       % x-position of receivers
    receiverz = 0.05 : .1 : 3;              % z-positions of receivers
    nsource = length(sourcez); nreceiver = length(receiverz);

    % Calculate acquisition geometry (multiple-offset gather)
    for j = 1 : nsource
        for i = 1 : nreceiver
            data( ( j - 1 ) * nreceiver + i , :) = [sourcex sourcez(j) receiverx receiverz(i)];
        end
    end

    % Calculate forward modeling kernel, Courtesy Dr. James Irving, UNIL
    Extra.J = tomokernel_straight(data,x,z); % Distance of ray-segment in each cell for each ray

    % Grid-cells in horizontal and vertical direction
    Extra.dimhor = length(x) - 1; Extra.dimver = length(z) - 1;

    Extra.por = 0.36;                       % Porosity field
    nz = 30;                                % Original discretization of water saturation model
    nx = 40;
    wcon = load('Sw.dat');                  % Load water saturation model: Courtesy of Dr. Michael Kowalsky (LBNL), Kowalsky et al. (2005; WRR)
    count = 0;
    for k = 1 : nz
        for i = 6 : nx - 5                  % Make model 3 by 3 m
            wcont( k , i - 5) = wcon( nz - k + 1 , i);
        end
    end
    wcont = wcont * Extra.por;              % Transform saturation into water content
    wcont = dct2(wcont);                    % Transform water content model into DCT space

    wtrunc = zeros(Extra.dimver,Extra.dimhor);
    for j = 1 : Extra.parz
        for i = 1 : Extra.parx
            wtrunc(j,i) = wcont(j,i);       % Truncate at same order as inverse parameterization
        end
    end
    wtrunc = idct2(wtrunc);                 % Do inverse DCT
    wtrunc = wtrunc';

    % Translate water content model into slowness using the Refractive Index/CRIM model
    Extra.pw = 81;                          % Permittivity of water
    Extra.pa = 1;                           % Permittivity of air
    Extra.ps = 5;                           % Permittivity of mineral grains
    slowtrue = wtrunc(:) * sqrt(Extra.pw) + (Extra.por-wtrunc(:)) * sqrt(Extra.pa) + (1 - Extra.por) * sqrt(Extra.ps); % sqrt of effective permittivity
    slowtrue = slowtrue / 0.3; % True slowness field

    % Simulate data for true slowness model
    % Add Gaussian uncorrelated noise with a standard deviation of 1 ns.
    Extra.datasim = Extra.J * slowtrue + Extra.error * randn(nsource * nreceiver,1); % Simulate data
    Extra.model_DCT = zeros(Extra.dimver,Extra.dimhor); % Matrix in which proposed model is assigned

    % Scale DCT coefficients such that all models are possible
    dum = zeros(Extra.dimver,Extra.dimhor);
    count = 1;
    for i = 1 : Extra.parz
        for j = 1 : Extra.parx
            if (i>1 | j>1)
                count = count + 1;
                dum(i,j) = 1;
                dummy = idct2(dum);
                ParRange.minn(count) = -1.7 / max(max(abs(dummy)));  % 0.2657
                ParRange.maxn(count) =  1.7 / max(max(abs(dummy)));  % 0.2657
                dum(i,j) = 0;
            end
        end
    end
    % ---------------------------------------------------------------------

    % Provide information to do initial sampling ("COV") --> The initial
    % chain positions are concentrated to the middle of the parameter ranges
    % This will speed up the calculation -- but cannot be done in practice!

    % Define mu and covariance of multinormal distribution used with "COV"
    MCMCPar.mu = ParRange.minn + 0.5 * ( ParRange.maxn - ParRange.minn );   % Provide mean of initial sample
    MCMCPar.cov = 0.001 * diag( ParRange.maxn - ParRange.minn );            % Initial covariance

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange);

end;

if example == 15,    % HYMOD rainfall - runoff model with measurement error estimation
    
    % ---------------------------- Check the following 3 papers ------------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of  %
    %      input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain %
    %      Monte Carlo simulation, Water Resources Research, 44, W00B09,                        %
    %      doi:10.1029/2007WR006720, 2008.                                                      %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson, Equifinality of formal    %
    %       (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling?, Stochastic %
    %       Environmental Research and Risk Assessment, 23(7), 1011-1026, 				        %
    %       doi:10.1007/s00477-008-0274-y, 2009                                                 %
    %                                                                                           %
    %   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution      %
    %       Metropolis algorithm for optimization and uncertainty assessment of hydrologic      %
    %       model parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003. %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.ndraw = 10000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 2;                        % Define likelihood function -- Sum of Squared Error

    % Define Func_name
    Func_name = 'hymodMATLAB';

    % Give the parameter ranges (minimum and maximum values) of the model 
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10 ]; 
    ParRange.maxn = [500 2.00 0.99 0.10 0.99 ];

    % ---------------------------------------------------------------------
    % Which measurement error shall we take? (with MCMCPar.lik = 2 or 7 !!)
    homoscedastic = 'No';
    
    % We attempt to estimate the Measurement error (homoscedastic)
    if strcmp(lower(homoscedastic),'yes'),
        % Sigma independent of magnitude of y (Measurement.MeasData)
        Measurement.Sigma = inline('a'); 
        % How many parameters does this error model have?
        Measurement.n = 1;
        % Now add "a" to the parameter ranges ("a" will be estimated)
        ParRange.minn = [ParRange.minn 0.5]; ParRange.maxn = [ParRange.maxn 100];
    elseif strcmp(lower(homoscedastic),'no'), % heteroscedastic
        % Sigma linearly dependent on y (Measurement.MeasData)
        Measurement.Sigma = inline('a * y + b'); 
        % "a" is slope, and "b" is intercept !! (make sure that Sigma > 0 )

        % How many parameters does this error model have?
        Measurement.n = 2;
        % Add "a" and "b" to parameter ranges (in order of "b" and "a" !!)
        ParRange.minn = [ParRange.minn 0 0]; ParRange.maxn = [ParRange.maxn 1 1];
    end;
    % NOTE: ONE CAN SPECIFY ANY TYPE OF HETEROSCEDASTIC ERROR MODEL! 
    % ---------------------------------------------------------------------
    
    % How many parameters do we have total?
    MCMCPar.n = size(ParRange.minn,2); 
    
    % Load the Leaf River data
    load bound.txt;

    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795;

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.E = bound(1:Extra.MaxT,5); Extra.P = sum(bound(1:Extra.MaxT,6:9),2);

    % Area factor to translate HYMOD output in mm/d to m3/s (calibration data); (area Leaf River is 1944 km2)
    Extra.F = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);

    % Define the measured streamflow data
    Measurement.MeasData = bound(65:Extra.MaxT,4);

    % We need to specify the Measurement error of the data in Measurement.Sigma
    % With option 3, Measurement.Sigma is integrated out the likelihoon function
    % With any other option, Sigma needs to be defined

    % We can estimate the measurement error directly if we use temporal differencing
    % The function MeasError provides an estimate of error versus flow level
    % out = MeasError(Measurement.MeasData;
    % For the Leaf River watershed this results in a heteroscedastic error
    % that is about 10% of the actual measured discharge value, thus
    % You can check this by plotting out(:,1) versus out(:,2)
    % Measurement.Sigma = 0.1*Measurement.MeasData; % --> option 2

    % Run the MT-DREAM_ZS algorithm
    [chain,X,Z,output,fx] = mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);

end;

% Create a single matrix with values sampled by chains
ParSet = GenParSet(chain,MCMCPar);

% --------------------------------------------------------------------------------------------- %
% ------------------------------------ POSTPROCESSING ----------------------------------------- %
%                                                                                               %
% For postprocessing of results --> please go to directory \PostProcessing and run the          %
% "postprocMCMC" script. This will compute various statistic and create a number of different   %
% plots, including the R_stat convergence diagnostic, marginal posterior parameter              %
% distributions, two-dimensional correlation plots of the posterior parameter samples, and      %
% parameter and total uncertainty posterior simulation uncertainty ranges.                      %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %
% --------------------------------------------------------------------------------------------- %