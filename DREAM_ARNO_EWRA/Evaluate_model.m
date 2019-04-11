function [fx] = Evaluate_model(x,MCMCPar,Measurement,Func_name,Extra);
% This function computes the likelihood and log-likelihood of each d-vector
% of x values
%
% Code both for sequential and parallel evaluation of model ( = pdf )
%
% Written by Jasper A. Vrugt

global EXAMPLE_dir;

% Check whether to store the output of each model evaluation (function call)
if ( strcmp(lower(MCMCPar.modout),'yes') ) && ( Measurement.N > 0 ),
    
    % Create initial fx of size model output by MCMCPar.seq
    fx = NaN(Measurement.N,MCMCPar.seq);
    
end;

% Now create the function handle (transparancy needed and thus in function with parfor loop!!!)
f_handle = eval(['@(x)',char(Func_name),'(x,Extra)']);

% Now go to example directory
% cd(EXAMPLE_dir)

% Now evaluate the model
if ( MCMCPar.CPU == 1 )         % Sequential evaluation
    
    % Loop over each d-vector of parameter values of x using 1 worker
    for ii = 1:size(x,1),
        
        % Execute the model and return the model simulation
        fx(:,ii) = f_handle(x(ii,:));
        
    end;
    
elseif ( MCMCPar.CPU > 1 )      % Parallel evaluation
    
    % If IO writing with model --> worker needs to go to own directory
    if strcmp(MCMCPar.IO,'Yes'),
        
        % Loop over each d-vector of parameter values of x using N workers
        parfor ii = 1:size(x,1),
            
            % Determine work ID
            t = getCurrentTask();
                        
              % Go to right directory (t.id is directory number)
            evalstr = strcat(EXAMPLE_dir,'\',num2str(t.id)); cd(evalstr)
            
            % Execute the model and return the model simulation
            fx(:,ii) = f_handle(x(ii,:));
            
        end;
        
    else
        
        % Loop over each d-vector of parameter values of x using N workers
        parfor ii = 1:size(x,1),
            
            % Execute the model and return the model simulation
            fx(:,ii) = f_handle(x(ii,:));
            
        end;
        
    end;
    
end;

% Now make sure we are back in DREAM directory
% cd(DREAM_dir);