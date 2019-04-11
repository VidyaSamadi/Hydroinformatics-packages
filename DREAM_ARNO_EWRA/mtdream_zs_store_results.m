function [ T ] = mtdream_zs_store_results ( MCMCPar , fx , Measurement , id );
% Stores the results of DREAM to binary files

% If T generations have been done, store information

%     % Append current states of X to chain files "chain_xx.bin"
%     for ii = 1 : MCMCPar.seq,
%         % Open ii'th chain
%         evalstr = strcat('fid_',num2str(ii),' = fopen(''chain_',num2str(ii),'.bin'',''',num2str(id),''',''n'');'); eval(evalstr);
%         % Now write X(ii,:) to ii'th chain
%         evalstr = strcat('fwrite(fid_',num2str(ii),',','X(ii,1:MCMCPar.n+2)''',',','''double''',');'); eval(evalstr);
%         % Now close file
%         evalstr = strcat('fclose(fid_',num2str(ii),');'); eval(evalstr);
%     end;

% Append current model simulations of X to file "fx.bin"
if ( strcmp(lower(MCMCPar.modout),'yes') ) && ( Measurement.N > 0 )
    % Now open the file to append new simulations
    evalstr = strcat('fid_fx = fopen(''fx.bin'',''',num2str(id),''',''n'');'); eval(evalstr);
    % Now append
    fwrite(fid_fx,fx,'double');
    % Now close file
    fclose(fid_fx);
end;

% And set the counter to 0
T = 0;