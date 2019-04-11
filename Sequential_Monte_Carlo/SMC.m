%% Sequential Monte Carlo
% Definitions:
    % Model used: Australian Water Balance Model (AWBM)
    % nter = number of particles (N)
    % n = length of data with initial burn in 
    % data_length = length of data without initial burn in. (Bates and Campbell, 2001)
    % dimen = number of parameters
    % ndimen = column for the parameter variance
    % ESS = Effective Sample Size
    % gammax_ss = dynamic sequential distribution
    % A1, A2, A3 = AWBM spatial areas
    % C1, C2, C3 = AWBM surface storages
    % K = daily recession constant
    % BFI = base flow index
    % V = parameter variance
    
% Inputs = prt:Percipitation (daily), et:Evapotranspiration (daily), ot:Observed Runoff (daily)
% Output = matrix of parameters (nter x dimen)  
    % The columns are arranged accordingly as below;
        % C1 C2 C3 A1 A2 K V BFI
    % ESS values per iteration 
    % SMC weights
    % maximum a posteri
% Catchment's information:
    % Never Never River, a 51km2 catchment located at the Gleniffer Bridge
    % in New South Wales, Australia.

    
    

    
%% Glabal variables (will be used repeteadly in the functions)
global nter n data_length dimen ndimen prt et ot

%% Load the  catchment data
[h,d]=hdrload('205014modrun.dat');
data=d;
prt=data(1:5151,4);
ot=data(1:5151,5);
et=data(1:5151,6);
ot=ot(151:length(prt)); % removing the first 150 days of data

%% Definition settings

nter=1000; % used to be 1000 particles
n=length(prt);
data_length=length(prt)-150;
dimen=8;
ndimen=dimen-1;

%% Given the prior mean and variance
[beta_mean,beta_var]=modelA(1.4,2.6);
A1m=beta_mean;
A1v=beta_var;
[beta_mean,beta_var]=modelA(2,2.5);
A2m=beta_mean;
A2v=beta_var;
[mbeta_mean,mbeta_var]=modelMB(0.41,1.5,5.19,9.74,8.24);
BFIm=mbeta_mean;
BFIv=mbeta_var;
[mbeta
  [prior]=priorWeibull(2.16,204,S3m);
  p1C3=prior;
  [prior]=priorMixBeta(0.271,51.9,4.17,255,9.6,Km);
  p1K=prior;
  [prior]=priorChi(8.6,46,Vm);
  p1var=prior;
  [prior]=priorMixBeta(0.41,1.5,5.19,9.74,8.24,BFIm);
  pbfi=prior;
  pdfx(jj)=log(p1C1)+log(p1C2)+log(p1C3)+log(p1A1)+log(p1A2)+log(p1K)...
           +(p1var)+log(pbfi);
  
end

%% Initialize the SMC sampler
S=1; % start of the distributional sequence
Send=1000;
PP=3.11; % exponential factor to monitor the length of the distributional sequence
Sx=1; % monitor the number of iteartions in the distributional path 
Sg=1; % A higher Sg, creates a smaller PP: reducing the change between each sequence
ESS_mat=zeros(Send,1);
w=repmat((1/nter),1,nter);
pre=1;
gammax_ss=1/10000000; % initial value of the distributional sequence.
gammax_ss_old=0;
g_ss=zeros(Send,1);
MAP=zeros(Send,1);
WNSE=zeros(Send,1);

if gammax_ss <= 1
  while gammax_ss <=1
        
    %% Generating proposed parameters
    den1_new=zeros(nter,dimen);
    for i = 1:nter
       
        den1_new(i,:)=mvnrnd(den1(i,:),cov_tita1);
        if  den1_new(i,1)<=0 || den1_new(i,2)<=0 || den1_new(i,3)<=0 || den1_new(i,4)+den1_new(i,5)>=1 ||...
            den1_new(i,4)<=0 || den1_new(i,5)<=0 || den1_new(i,6)>=1 || den1_new(i,6)<=0 || den1_new(i,7)<=0 || den1_new(i,8)<=0 || den1_new(i,8)>=1 
        
            while den1_new(i,1)<=0 || den1_new(i,2)<=0 || den1_new(i,3)<=0 || den1_new(i,4)+den1_new(i,5)>=1 ||...
                  den1_new(i,4)<=0 || den1_new(i,5)<=0 || den1_new(i,6)>=1 || den1_new(i,6)<=0 || den1_new(i,7)<=0 || den1_new(i,8)<=0 || den1_new(i,8)>=1 
         
              den1_new(i,:)=mvnrnd(den1(i,:),cov_tita1);
            end
        end
    end
    
    %% MCMC step
    den1_newx=zeros(nter,dimen);
    Pri_new=zeros(nter,1);

    for j= 1:nter
    
    [qt]=AWBM(den1(j,1),den1(j,2),den1(j,3),den1(j,4),den1(j,5),den1(j,6),den1(j,8));
     Q1=qt(151:length(prt));
            
    [qt]=AWBM(den1_new(j,1),den1_new(j,2),den1_new(j,3),den1_new(j,4),den1_new(j,5),den1_new(j,6),den1_new(j,8));
     Q1_new=qt(151:length(prt));
     
     % prior values
     [prior]=priorWeibull(2.16,68,den1(j,1));
     pC1=prior;
     [prior]=priorWeibull(2.16,68,den1_new(j,1));
     p0C1=prior;
     [prior]=priorWeibull(2.16,102,den1(j,2));
     pC2=prior;
     [prior]=priorWeibull(2.16,102,den1_new(j,2));
     p0C2=prior;
     [prior]=priorWeibull(2.16,204,den1(j,3));
     pC3=prior;
     [prior]=priorWeibull(2.16,204,den1_new(j,3));
     p0C3=prior;
     
     [prior]=priorBeta(1.4,2.6,den1(j,4));
     pA1=prior;
     [prior]=priorBeta(1.4,2.6,den1_new(j,4));
     p0A1=prior;
     [prior]=priorBeta(2,2.5,den1(j,5));
     pA2=prior;
     [prior]=priorBeta(2,2.5,den1_new(j,5));
     p0A2=prior;

     [prior]=priorMixBeta(0.271,51.9,4.17,255,9.6,den1(j,6));
     pK=prior;
     [prior]=priorMixBeta(0.271,51.9,4.17,255,9.6,den1_new(j,6));
     p0K=prior;
     
     [prior]=priorChi(8.6,46,den1(j,7));
     pv=prior;
     [prior]=priorChi(8.6,46,den1_new(j,7));
     p0v=prior;
     
     [prior]=priorMixBeta(0.41,1.5,5.19,9.74,8.24,den1(j,8));
     pBFI=prior;
     [prior]=priorMixBeta(0.41,1.5,5.19,9.74,8.24,den1_new(j,8));
     p0BFI=prior;
     
          
    [den1_accept,prior_accept]=modelacpsmc(Q1,Q1_new,den1(j,:),den1_new(j,:),pC1,pC2,pC3,pA1,pA2,pK,pv,pBFI,p0C1,p0C2,p0C3,p0A1,p0A2,p0K,p0v,p0BFI);
    den1_newx(j,:)=den1_accept;
    Pri_new(j)=prior_accept;
   
   
    end
    den1_acc=den1_newx;

    %% Weighting
    weight=zeros(nter,1);
    cal_w=zeros(nter,1);
    for m = 1:nter

        [qt]=AWBM(den1_acc(m,1),den1_acc(m,2),den1_acc(m,3),den1_acc(m,4),den1_acc(m,5),den1_acc(m,6),den1_acc(m,8));
             Q1_new=qt(151:length(prt));

        vv=den1_acc(m,ndimen);     
        alp =  gammax_ss*((Pri_new(m))+((-(data_length/2))*log(2*pi*vv) + sum(-((ot-Q1_new).^2)/(2*vv))));
        prix = (1-gammax_ss)*pdfx(m);
        alp_old =  gammax_ss_old*((Pri_new(m))+((-(data_length/2))*log(2*pi*vv) + sum(-((ot-Q1_new).^2)/(2*vv))));
        prix_old = (1-gammax_ss_old)*pdfx(m);
        weight(m) = alp + prix;
        weight_old = alp_old + prix_old;
        cal_w(m) = weight(m) - weight_old;

    end

    w= w .* exp(cal_w)';  
    w= w ./ sum(w);
    % Calculating the ESS
    ESS=(sum(w))^2/(sum(w.^2));
    ESS_mat(S)=ESS;
    g_ss(S)=gammax_ss;
    gammax_ss_old=gammax_ss;
    

    %% NSE
    nse=zeros(nter,1);
    like=zeros(nter,1);
    for m2 = 1:nter    
       [qt]=AWBM(den1_acc(m2,1),den1_acc(m2,2),den1_acc(m2,3),den1_acc(m2,4),den1_acc(m2,5),den1_acc(m2,6),den1_acc(m2,8));
         Q=qt(151:length(prt));
        
        vv=den1_acc(m2,ndimen);
        nse(m2) = w(m2)*(1-((sum((ot-Q).^2))/(sum((ot-mean(ot)).^2))));
        like(m2) = w(m2)*((Pri_new(m2))+((-(data_length/2))*log(2*pi*vv) + sum(-((ot-Q).^2)/(2*vv))));


    end
    WNSE(S) = sum(nse);
    %% Resampling
    %% Dynamic Distributional Sequence
    %% Threshold if the sampler collpase to NaN
    if isnan(ESS)
    
       Pcount=0.1488*exp(-0.2127*Sg); % Recalculate PP
       PP=PP-Pcount;
       gammax_ss=(7*10^-8)*((Sx)^PP); 
       gammax_ss_old=(7*10^-8)*((Sx-1)^PP);
       Sg=Sg+1;

       ind=randsample(1:nter,nter,true,w_old);
       den1=den1_old(ind,:);
       w=repmat((1/nter),1,nter);

       cov_tita1=covweight(den1,w);

       Sx=Sx+1;
       pre=pre+1;
       S=S+1;
   
    else
        %% Threshold for resampling
        if ESS< nter/2 

            ind=randsample(1:nter,nter,true,w);
            den1=den1_acc(ind,:);
            w=repmat((1/nter),1,nter);
            %% Second threshold if the sampler collapse        
            if ESS < 0.4*(nter/2)
                Pcount=0.1488*exp(-0.2127*Sg);
                PP=PP-Pcount;
                gammax_ss_old=(7*10^-8)*((Sx-1)^PP);
                Sg=Sg+1;

                ind=randsample(1:nter,nter,true,w_old);
                den1=den1_old(ind,:);
                w=repmat((1/nter),1,nter);
            end

           cov_tita1=covweight(den1,w);
           gammax_ss=(7*10^-8)*((Sx)^PP);

        else
            %% If no resampling
            den1 = den1_acc;
            cov_tita1=covweight(den1,w);
            gammax_ss=(7*10^-8)*((Sx)^PP);

        end

            w_old=w;
            den1_old=den1;

            Sx=Sx+1;
            pre=pre+1;
            S=S+1;
    
    end

  end
end

S = S-1
MAP = sum(like)
mean(den1)
CI95 = quantile(den1,[.025 .975])


xx=[0 S];
yy=[nter/2 nter/2];
plot(ESS_mat(1:S))
line(xx,yy)
subplot(2,1,1)
plot(WNSE(1:S))
subplot(2,1,2)
plot(ESS_mat(1:S))
line(xx,yy)



