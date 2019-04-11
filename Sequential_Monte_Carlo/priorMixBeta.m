%% Prior : Mix-Beta distribution
function[varargout]=priorMixBeta(tao,a1,b1,a2,b2,tita)

    p_mix1 = (gammaln(a1+b1)-(gammaln(a1)+gammaln(b1)))+log(tita^(a1-1))+log((1-tita)^(b1-1));
    p_mix2 = (gammaln(a2+b2)-(gammaln(a2)+gammaln(b2)))+log(tita^(a2-1))+log((1-tita)^(b2-1));
    prior = (tao*(exp(p_mix1)))+((1-tao)*(exp(p_mix2)));
    
varargout{1}=prior;


