function [MCMCPar] = mtdream_zs_calc_setup(MCMCPar);
% Sets up sequential / parallel

global DREAM_dir EXAMPLE_dir;

% Now check if we want parallel execution of chains or not?
if strcmp( lower ( MCMCPar.parallel ) , 'no' )
    
    % We use 1 CPU (processor)
    MCMCPar.CPU = 1;

else
    
    % First close possible parallel jobs
    try, matlabpool close force local; catch dummy = 0, end

    % Now open available cores
    matlabpool open
    
    % Check whether computer has multiple cores?
    MCMCPar.CPU = matlabpool('size'); 
    
    % If input/output writing is done we need directories for each worker
    if strcmp(MCMCPar.IO,'Yes'),
        
        % Go to directory with problem files
        cd(EXAMPLE_dir)
        
        % Create first part of copy expression
        a = strcat('copy "'); 
        
        % Create the directories
        for ii = 1:min(MCMCPar.CPU,MCMCPar.mt * MCMCPar.seq),
            
            % Create the directories
            mkdir(strcat(num2str(ii)));
            
            % And copy the files
            b = strcat(a,EXAMPLE_dir,'\*.*"',{' '},'"',EXAMPLE_dir,'\',strcat(num2str(ii)),'"'); dos(char(b));
            
        end;
        
        cd(DREAM_dir);
    
    end;

end;