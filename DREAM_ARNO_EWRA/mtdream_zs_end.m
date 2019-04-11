function [chain,Z,output,fx] = mtdream_zs_end(MCMCPar,Measurement,chain,output,iteration,iloc,pCR,fid_Z,argout);
% Variables have been pre-allocated --> need to remove zeros at end

global DREAM_dir EXAMPLE_dir

% Start with CR
output.CR = output.CR(1:iteration-1,1:size(pCR,2)+1);
% Then R_stat
output.R_stat = output.R_stat(1:iteration-1,1:MCMCPar.n+1);
% Then AR
output.AR = output.AR(1:iteration-1,1:2);
% Adjust last value (due to possible sudden end of for loop)
% output.AR(counter,1:2) = [output.AR(counter-1,1) + MCMCPar.seq  mean(output.AR(2:counter-1,2))];
% Then chain
chain = chain(1:iloc,1:MCMCPar.n+2,1:MCMCPar.seq);
% Then the history Z
%Z = Z(1:MCMCPar.m,1:MCMCPar.n+2);

% First close the file
fclose(fid_Z);
% Open the binary file with thinned sample history
fid_Z = fopen('Z.bin','r','n');
% Now read the binary file
Z = fread(fid_Z, [MCMCPar.n+2,MCMCPar.m],'double')';
% Now close the file again
fclose(fid_Z);

% Now calculate the convergence diagnostics for individual chains using the CODA toolbox
fid = fopen('MT-DREAM_ZS_diagnostics.txt','w');
% Now loop over each chain
for j = 1:MCMCPar.seq,
    % First calculate diagnostics
    diagnostic{j} = coda(chain(floor(0.5*iloc):iloc,1:MCMCPar.n,j));
    % Now write to file DREAM.out
    diagnostic{j}(1).chain_number = j; prt_coda(diagnostic{j},[],fid);
end;
% Now close the file again
fclose(fid);

% If five output arguments are requested then return fx
if argout == 5,
    switch lower(MCMCPar.modout)
        case 'no'
            % Return an empty matrix
            fx = [];
        case 'yes'
            % Open the binary file with model simulations
            fid_fx = fopen('fx.bin','r','n');
            % Now read the binary file
            fx = fread(fid_fx, [Measurement.N,MCMCPar.ndraw],'double')';
            % Now close the file again
            fclose(fid_fx);
    end;
end;

% Close MATLAB pool (if CPU > 1) and remove file if
if MCMCPar.CPU > 1,
    % Close the matlab pool
    matlabpool('close');
    % If input output writing, then remove directories
    if strcmp(lower(MCMCPar.IO),'yes');
        % Go to directory with problem files
        cd(EXAMPLE_dir)
        % Remove the directories
        for ii = 1:min(MCMCPar.CPU,MCMCPar.seq)
            % Remove each worker directory
            rmdir(strcat(num2str(ii)),'s');
        end;
        cd(DREAM_dir)
    end;
end;

% Remove path of dream
rmpath(DREAM_dir);