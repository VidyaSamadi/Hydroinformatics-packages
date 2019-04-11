function s = LHSU(xmin,xmax,nsample)
% Generates sample using Latin Hypercube Sampling
% Adapted from Budiman (2003)

% Define the number of variables
nvar = length(xmin);

% Initialize with random elements
ran = rand(nsample,nvar);

% Initialize output array with zeros
s = zeros(nsample,nvar);

% Loop over number of variables
for j = 1:nvar
   idx = randperm(nsample);
   P = (idx'-ran(:,j))/nsample;
   s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
end