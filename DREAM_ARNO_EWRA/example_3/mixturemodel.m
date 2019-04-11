function [p] = mixturemodel(x,Extra);
% Example of mixture model with two multinormal distributions centered at -5, 5 and 15. 

% Define the sigma of the target distribution
sigma = eye(size(x,2));

% Determine scaling factor
F = 1e12;

% Modes at -5, 5 and 15

% Calculate probability density at given value of x
% p = F^(size(x,2)) * (1/6 * mvnpdf(x, -5 , sigma ) + 2/6 * mvnpdf(x, 5 , sigma ) + 3/6 * mvnpdf(x, 15 , sigma ))

% Solve differently to avoid numerical problems! (zero probability)
p1 = F * normpdf( x, -5, 1); p2 = F * normpdf( x,  5, 1); p3 = F * normpdf( x, 15, 1);

% Compute the sum of the mixtures
p = max ( sum( 1/6 * prod(p1) + 2/6 * prod(p2) + 3/6 * prod(p3) ) , 1e-299 );