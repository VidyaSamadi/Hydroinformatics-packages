function [MCMCPar,pCR,CR,lCR,Iter,counter,teller,new_teller,hist_logp,Sequences,Table_JumpRate,...
    Reduced_Seq,iloc,output] = InitVariables(MCMCPar,Extra)
% Initializes important variables for use in the algorithm

% Calculate the parameters in the exponential power density function of Box and Tiao (1973)
[MCMCPar.Cb,MCMCPar.Wb] = CalcCbWb(MCMCPar.Gamma); 

% Initialize the array that contains the history of the log_density of each chain
hist_logp = zeros(floor(MCMCPar.ndraw/MCMCPar.seq),MCMCPar.seq);

% Derive the number of elements in the output file
Nelem = floor(MCMCPar.ndraw/MCMCPar.seq) + 1;

% Initialize output information -- AR
output.AR = zeros(Nelem,2); 
output.AR(1,1:2) = [MCMCPar.seq -1];

% Initialize output information -- Outlier chains
output.outlier = [];

% Initialize output information -- R statistic
output.R_stat = zeros(floor(Nelem/MCMCPar.steps),MCMCPar.n+1);

if strcmp(Extra.pCR,'Update'),
    % Calculate multinomial probabilities of each of the nCR CR values
    pCR = (1/MCMCPar.nCR) * ones(1,MCMCPar.nCR);
    % Calculate the actual CR values based on p
    [CR] = GenCR(MCMCPar,pCR); lCR = zeros(1,MCMCPar.nCR);
else
    pCR = 1/MCMCPar.nCR; lCR = [];
    % Define 
    [CR] = pCR * ones(MCMCPar.seq,MCMCPar.steps); lCR = MCMCPar.seq * MCMCPar.steps;
end;

% Initialize output information -- N_CR  
output.CR = zeros(floor(Nelem/MCMCPar.steps),size(pCR,2)+1);

if strcmp(Extra.save_in_memory,'Yes');
    % Initialize Sequences with zeros
    Sequences = zeros(floor(1.25 * Nelem),MCMCPar.n+1,MCMCPar.seq);
else
    Sequences = [];
end;

% Generate the Table with JumpRates (dependent on number of dimensions and number of pairs)
for zz = 1:MCMCPar.DEpairs,
    Table_JumpRate(:,zz) = 2.38./sqrt(2 * zz * [1:MCMCPar.n]'); 
end;

% Check whether will save a reduced sample
if strcmp(Extra.reduced_sample_collection,'Yes');
    % Initialize Sequences with zeros
    Reduced_Seq = zeros(floor(Nelem/Extra.T),MCMCPar.n+2,MCMCPar.seq);
else
    Reduced_Seq = [];
end;

% Initialize Iter and counter
Iter = MCMCPar.seq; counter = 2; iloc = 1; teller = 2; new_teller = 1; 

% Change MCMCPar.steps to make sure to get nice iteration numbers in first loop
MCMCPar.steps = MCMCPar.steps - 1;