%% Bayesian Inference
function[varargout] = modelacpsmc(Q1,Q1_new,den1,den1_new,pC1,pC2,pC3,pA1,pA2,pK,pv,pBFI,p0C1,p0C2,p0C3,p0A1,p0A2,p0K,p0v,p0BFI)
%% Identify the Global variables
global data_length ot ndimen

%% Parameter variances
   vv1=den1(ndimen);
   vv1_new=den1_new(ndimen);

%% Joint prior distributions  
   prior_old=log(pC1)+log(pC2)+log(pC3)+log(pA1)+log(pA2)+log(pK)+(pv)+log(pBFI);
   
   prior_new=log(p0C1)+log(p0C2)+log(p0C3)+log(p0A1)+log(p0A2)+log(p0K)+(p0v)+log(p0BFI);
   
   prior_old=real(prior_old);
   prior_new=real(prior_new);
%% Likelihoods   
   likelihood_old=(-(data_length/2)*log(2*pi*vv1)+sum(-((ot-Q1).^2)/(2*vv1)));
   likelihood_new=(-(data_length/2)*log(2*pi*vv1_new)+sum(-((ot-Q1_new).^2)/(2*vv1_new)));
   
%% Acceptance probabilities   
alpha=exp(likelihood_new + prior_new - likelihood_old - prior_old);

if isnan(alpha)
    alpha=0;
end

U=unifrnd(0,1,1);

if alpha>U 
    den1_accept=den1_new;
    prior_accept=prior_new;   
else 
    den1_accept=den1;
    prior_accept=prior_old;
    
end    

varargout{1} = den1_accept;
varargout{2} = prior_accept;

