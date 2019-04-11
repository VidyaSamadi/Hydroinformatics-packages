function [stop,fid] = mtdream_zs_check(MCMCPar,Measurement);
% Check for setup errors

% Assign stop to be No
stop = 'No';

% First close all files
fclose('all'); 

% Then delete the diagnostics output file
dos('del MT-DREAM_ZS_diagnostics.txt');

% open an output file with warnings
fid = fopen('warning_file.txt','w');
fprintf(fid,'-------------- DREAM_{(ZS)} warning file --------------\n');

% Check number of chains
if MCMCPar.seq < 3,
    % Warning -- not enough chains to do sampling -- increase number of chains!
    evalstr = char(strcat('DREAM WARNING: Inufficient number of chains -> Use at least MCMCPar.seq = ',{' '},num2str((2 * MCMCPar.DEpairs) + 1),{' '},'chains \n'));   
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
    % Stop DREAM
    stop = 'Yes';
end;

% Check that MCMCPar.m is large enough to create offspring
if MCMCPar.m0 < (2 * MCMCPar.DEpairs * MCMCPar.seq),
    % Warning -- not enough chains to do sampling -- increase number of chains!
    evalstr = char(strcat('DREAM_{(ZS)} WARNING: size of Z not sufficient -> Increase MCMCPar.m0 to at least ',{' '},...
        num2str(2 * MCMCPar.DEpairs * MCMCPar.seq),'\n'));   
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
end;

% Check whether we specified the measurement Sigma
if ( MCMCPar.lik == 2 || MCMCPar.lik == 7 ) && ( isfield(Measurement,'Sigma') == 0 ),
    % Warning -- Measurement.Sigma needs to be specified!!
    evalstr = char('DREAM WARNING: Measurement.Sigma needs to be specified either as inline function or one or multiple numerical values!!\n');   
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
    % Stop DREAM
    stop = 'Yes';
end;

% Check whether the length of the user specified measurement Sigma is correct
if ( isfield(Measurement,'Sigma') == 1 ) && ( isfield(Measurement,'n') == 0 ),
    if (prod(size(Measurement.Sigma)) ~= prod(size(Measurement.MeasData))) && ( prod(size(Measurement.Sigma)) > 1 )
        % Warning -- Measurement.Sigma incorrect length!!
        evalstr = char('DREAM WARNING: Heteroscedastic error, but length of Measurement.Sigma is not equal to that of the observations!!\n');   
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Stop DREAM
        stop = 'Yes';
    elseif ( sum(Measurement.Sigma<=0) > 0 ),
        % Warning -- Measurement.Sigma is negative!!
        evalstr = char('DREAM WARNING: At least one value of the specified Measurement.Sigma is negative or zero!!\n');   
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Stop DREAM
        stop = 'Yes';
    end;
end;