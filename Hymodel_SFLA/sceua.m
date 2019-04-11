%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bestx,bestf] = sceua(x0,bl,bu,maxn,kstop,pcento,peps,ngs,iseed,iniflg)
tic
% This is the subroutine implementing the SCE algorithm,
% written by Q.Duan, 9/2004
%
% taarif:
%  x0 = moteqayer;
%     = parametre behine shode dar enteha;
%  f0 = tabe hadaf
%     = maqadire tabe hadaf motabeqe parametr haye behine shode
%  bl = min meqdar motegyer
%  bu = max meqdar motegyer
%  iseed = random number (for repetetive testing purpose)
%  iniflg = alamat baraye moteqayer (=1, jozve jameiat avalie; dar qeire in soorat, shamel nist)
%  ngs = tedad zirmajmooeha
%  npg = tadad azae zirmajmooe
%  nps = number of members in a simplex
%  nspl = tedad jahesh haye arzyabi qabl az beham amikhtan
%  mings = min majmooe hay lazem da tey farayand behinesazi
%  maxn = max tedad tavabe arzyabi mojaz dar tey farayand
%  behinesazi
%  kstop = tedad tekrar halqe qabl az hamgerae
%  percento = darsad taqirat mojaz dar kstop loop qabl az hamgerae

% list motaqayer hay mahali
%    x(.,.) = mokhtasat noqat dar jameiat
%    xf(.) = maqadire tabe x(.,.)
%    xx(.) =mokhtasat yek noqte tanha x
%    cx(.,.) = mokhtasat noqat dar yek majmooe
%    cf(.) = maqadire tabe cx(.,.)
%    s(.,.) = mokhtasat noqat dar simplex konooni
%    sf(.) = function values of s(.,.)
%    bestx(.) = behtarin noqte da hamin halqe shuffling loop
%    bestf = function value of bestx(.)
%    worstx(.) = worst point at current shuffling loop
%    worstf = function value of worstx(.)
%    xnstd(.) = standard deviation of parameters in the population
%    gnrng = normalized geometri%mean of parameter ranges
%    lcs(.) = indices locating position of s(.,.) in x(.,.)
%    bound(.) = bound on ith variable being optimized
%    ngs1 = number of complexes in current population
%    ngs2 = number of complexes in last population
%    iseed1 = current random seed
%    criter(.) = vector containing the best criterion values of the last
%                10 shuffling loops

global BESTX BESTF ICALL PX PF Spar Qobs num sst minpar95 maxpar95 mintot95 maxtot95

% Initialize SCE parameters:
nopt=length(x0);
npg=2*nopt+1;
nps=nopt+1;
nspl=npg;
mings=ngs;
npt=npg*ngs;

bound = bu-bl;

% Create an initial population to fill array x(npt,nopt):
rand('seed',iseed)
x=zeros(npt,nopt);
for i=1:npt;
    x(i,:)=bl+rand(1,nopt).*bound;
end;

if iniflg==1; x(1,:)=x0; end;

nloop=0;
icall=0;
for i=1:npt;
    xf(i) = NS(nopt,x(i,:));
    icall = icall + 1;
end;
f0=xf(1);

% Sort the population in order of increasing function values;
[xf,idx]=sort(xf);
x=x(idx,:);

% Record the best and worst points;
bestx=x(1,:); bestf=xf(1);
worstx=x(npt,:); worstf=xf(npt);
BESTF=bestf; BESTX=bestx;ICALL=icall;

% Compute the standard deviation for each parameter
xnstd=std(x);

% Computes the normalized geometric range of the parameters
gnrng=exp(mean(log((max(x)-min(x))./bound)));

disp('The Initial Loop: 0');
disp(['BESTF  : ' num2str(bestf)]);
disp(['BESTX  : ' num2str(bestx)]);
disp(['WORSTF : ' num2str(worstf)]);
disp(['WORSTX : ' num2str(worstx)]);
disp(' ');

% Check for convergency;
if icall >= maxn;
    disp('*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT');
    disp('ON THE MAXIMUM NUMBER OF TRIALS ');
    disp(maxn);
    disp('HAS BEEN EXCEEDED.  SEARCH WAS STOPPED AT TRIAL NUMBER:');
    disp(icall);
    disp('OF THE INITIAL LOOP!');
end;

if gnrng < peps;
    disp('THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE');
end;

% Begin evolution loops:
nloop = 0;
criter=[];
criter_change=1e+5;

while icall<maxn & gnrng>peps & criter_change>pcento;
    nloop=nloop+1;

    % Loop on complexes (sub-populations);
    for igs = 1: ngs;

        % Partition the population into complexes (sub-populations);
        k1=1:npg;
        k2=(k1-1)*ngs+igs;
        cx(k1,:) = x(k2,:);
        cf(k1) = xf(k2);

        % Evolve sub-population igs for nspl steps:
        for loop=1:nspl;

            % Select simplex by sampling the complex according to a linear
            % probability distribution
            lcs(1) = 1;
            for k3=2:nps;
                for iter=1:1000;
                    lpos = 1 + floor(npg+0.5-sqrt((npg+0.5)^2 - npg*(npg+1)*rand));
                    idx=find(lcs(1:k3-1)==lpos); if isempty(idx); break; end;
                end;
                lcs(k3) = lpos;
            end;
            lcs=sort(lcs);

            % Construct the simplex:
            s = zeros(nps,nopt);
            s=cx(lcs,:); sf = cf(lcs);

            [snew,fnew,icall]=cceua(s,sf,bl,bu,icall,maxn);

            % Replace the worst point in Simplex with the new point:
            s(nps,:) = snew; sf(nps) = fnew;

            % Replace the simplex into the complex;
            cx(lcs,:) = s;
            cf(lcs) = sf;

            % Sort the complex;
            [cf,idx] = sort(cf); cx=cx(idx,:);

            % End of Inner Loop for Competitive Evolution of Simplexes
        end;

        % Replace the complex back into the population;
        x(k2,:) = cx(k1,:);
        xf(k2) = cf(k1);

        % End of Loop on Complex Evolution;
    end;

    % Shuffled the complexes;
    [xf,idx] = sort(xf); x=x(idx,:);
    PX=x; PF=xf;

    % Record the best and worst points;
    bestx=x(1,:); bestf=xf(1);
    worstx=x(npt,:); worstf=xf(npt);
    BESTX=[BESTX;bestx]; BESTF=[BESTF;bestf];ICALL=[ICALL;icall];

    % Compute the standard deviation for each parameter
    xnstd=std(x);

    % Computes the normalized geometric range of the parameters
    gnrng=exp(mean(log((max(x)-min(x))./bound)));

    disp(['Evolution Loop: ' num2str(nloop) '  - Trial - ' num2str(icall)]);
    disp(['BESTF  : ' num2str(bestf)]);
    disp(['BESTX  : ' num2str(bestx)]);
    disp(['WORSTF : ' num2str(worstf)]);
    disp(['WORSTX : ' num2str(worstx)]);
    disp(' ');

    % Check for convergency;
    if icall >= maxn;
        disp('*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT');
        disp(['ON THE MAXIMUM NUMBER OF TRIALS ' num2str(maxn) ' HAS BEEN EXCEEDED!']);
    end;

    if gnrng < peps;
        disp('THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE');
    end;

    criter=[criter;bestf];
    if (nloop >= kstop);
        criter_change=abs(criter(nloop)-criter(nloop-kstop+1))*100;
        criter_change=criter_change/mean(abs(criter(nloop-kstop+1:nloop)));
        if criter_change < pcento;
            disp(['THE BEST POINT HAS IMPROVED IN LAST ' num2str(kstop) ' LOOPS BY ', ...
                'LESS THAN THE THRESHOLD ' num2str(pcento) '%']);
            disp('CONVERGENCY HAS ACHIEVED BASED ON OBJECTIVE FUNCTION CRITERIA!!!')
        end;
    end;

    % End of the Outer Loops
end;

disp(['SEARCH WAS STOPPED AT TRIAL NUMBER: ' num2str(icall)]);
disp(['NORMALIZED GEOMETRIC RANGE = ' num2str(gnrng)]);
disp(['THE BEST POINT HAS IMPROVED IN LAST ' num2str(kstop) ' LOOPS BY ', ...
    num2str(criter_change) '%']);
Ntot = sst-1 
ModPred=SimRR(bestx);
% Compute the RMSE
RMSE = sqrt ( sum ( ( ModPred - Qobs).^2) / num);

% Now create the parameter and total uncertainty
for j = 1:Ntot
    Smod(:,j) = Spar(:,j) + normrnd(0,RMSE,num,1);
end;



% Now sort to get 95% ranges
for j = 1:num
    a = sort(Spar(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    minpar95(j,1) = a(floor(0.025 * Ntot)); maxpar95(j,1) = a(floor(0.975 * Ntot));
    % Same with total uncertainty
    a = sort(Smod(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    mintot95(j,1) = a(floor(0.025 * Ntot)); maxtot95(j,1) = a(floor(0.975 * Ntot));
end;

% Now plot the results 
plot([1:num],minpar95,[1:num],maxpar95); % These are the upper and lower parameter uncertainty bounds
hold on; plot([1:num],mintot95,'R',[1:num],maxtot95,'G'); % These are the upper and lower total uncertainty bounds
plot([1:num],Qobs,'r+');

% Now calculate percentage inside the 95% bound --> should be about 95%
Contained = 1 - length(find(Qobs<mintot95 |  Qobs> maxtot95))/num
Contained1 = 1 - length(find(Qobs<minpar95 |  Qobs> maxpar95))/num
toc
% Note the ranges do not look very good --> because of homoscedasticity. If
% the error is assumed to increase with flow level results look much better
return;
