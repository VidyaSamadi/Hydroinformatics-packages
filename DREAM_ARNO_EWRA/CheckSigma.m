function Measurement = CheckSigma(Measurement,MCMCPar);
% Now check how the measurement sigma is arranged (estimated or defined)

% First calculate the number of calibration data measurements
Measurement.N = size(Measurement.MeasData,1);

% Check whether the measurement error is estimated jointly with the parameters
if isfield(Measurement,'Sigma'),
    % Now check whether Sigma is estimated or whether user has defined Sigma
    if isfield(Measurement,'n'),
        % Set initial part of string
        str_sigma = strcat('Sigma = Measurement.Sigma('); 
        % Loop over n
        for j = 1:Measurement.n,
            % Define the input variables of the inline function in text
            evalstr = strcat('Measurement.a',num2str(j-1),'=','''x(ii,MCMCPar.n - ',num2str(j)-1,' )'';'); 
            % Now evaluate the text
            eval(evalstr); 
            % Now add to str_sigma
            str_sigma = strcat(str_sigma,'eval(Measurement.a',num2str(j-1),'),');
        end;
        % If Measurement.n > 1 --> Measured data is used
        if Measurement.n > 1,
            % Now finish the function call
            Measurement.str_sigma = strcat(str_sigma,'Measurement.MeasData);'); 
        else
            str_sigma = str_sigma(1:end-1); Measurement.str_sigma = strcat(str_sigma,');'); 
        end;
    end;
else
    % Sigma not estimated
end;