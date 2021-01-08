#
# R script for application of DREAM to SWAT applied to the Waccamaw watershed, the Carolinas.
#
#
#
#load dream. This line need to be executed if dream has already been loaded in the R session.
#
library(dream)
#
# Load and start the SNOW parallel processing library. This line need to be executed if snow
#has already been loaded in the R session. Also set the number of nodes equal to the number
#of MCMC chains to be used.
#
library(snow)
cl<-makeCluster(8,type="SOCK")
#
#the number of nodes and MCMC chains (nseq) is set at 8.
#The maximum number of simulations (ndraw) is limited to 24,000.
#The "snow.chains" method of parallelization should be used (parallel).
#The convergence criteria is left at its default value of control$Rthres=1.01.
#The modeler may consider limiting the number of simulations to a relatively low number until
#discovering how well the chains converge for the particular model and likelihood function.
#We found that when the Gaussian likelihood function was used for our swat model the 
#convergence criteria of control$R.thres=1.2 recommended by Vrugt et al (2009) was usually
#satisfied whitin 24,000 simulations,whereas with the laplace likelihood function it was not.
#However, best-stimates of parameters did not change significantly after 24,000
#simulations for either likelihood function/
#
control<-list(nseq=8,ndraw=40,parallel="snow.chains")
#
#The SWAT adjustment parameters and their ranges are shown below for our example. See pp.40-43
# of the SWAT-CUP4 user manual (Abbaspour, 2011) for additional guidance. These parameters must
#appear in the same order as they do in the "model.in.rows"statement.
#
pars<-list("r__CN2.mgt"=c(-0.50,0.50),"v__ALPHA_BF.gw"=c(0,1.0),"v__GW_DELAY.gw"=c(0,500),
"v__GW_REVAP.gw"=c(0.02,0.2),"v__RCHRG_DP.gw"=c(0,1),"r__SOL_K().sol"=c(-0.5,0.5),
"r__SOL_BD().sol"=c(-0.5,0.5),"v__ESCO.hru"=c(0.01,1),"v__EPCO.hru"=c(0.01,1),"v__SLSUBBSN.hru"=c(10,150),
"v__OV_N.hru"=c(0,0.8),"v__LAT_TTIME.hru"=c(0,180),"v__CH_K2.rte"=c(0,150),"v__CH_N2.rte"=c(0,0.3),
"r__SOL_AWC().sol"=c(-0.5,0.5),"r__GW_SPYLD.gw"=c(-0.5,0.5),"r__GWHT.gw"=c(-0.5,0.5),
"r__SHALLST.gw"=c(-0.5,0.5)
)
#
# Define the function that will be called by DREAM for every parameter set in every chain.
#
callSWAT<-function(id,pars){
model.in.rows<-c("r__CN2.mgt","v__ALPHA_BF.gw","v__GW_DELAY.gw","v__GW_REVAP.gw","v__RCHRG_DP.gw",
"r__SOL_K().sol","r__SOL_BD().sol","v__ESCO.hru","v__EPCO.hru","v__SLSUBBSN.hru","v__OV_N.hru",
"v__LAT_TTIME.hru","v__CH_K2.rte","v__CH_N2.rte","r__SOL_AWC().sol","r__GW_SPYLD.gw","r__GWHT.gw","r__SHALLST.gw")
model.in.file<-as.character(c("LREW1/model.in","LREW2/model.in",
			"LREW3/model.in","LREW4/model.in","LREW5/model.in",
			"LREW6/model.in","LREW7/model.in","LREW8/model.in"))
output.rch.file<-as.character(c("LREW1/output.rch","LREW2/output.rch",
			"LREW3/output.rch","LREW4/output.rch","LREW5/output.rch",
			"LREW6/output.rch","LREW7/output.rch","LREW8/output.rch"))
# Extract values from the observed data file to create the vector of observed values,
#and create vector of names of the batch files that will execute SWAT and SWAT_Edit.
 Q1o<-read.table("Qo.txt")[,1]
 Q2o<-read.table("Qo.txt")[,2]
 bfname<-as.character(c("runLREW1.bat","runLREW2.bat",
			"runLREW3.bat","runLREW4.bat","runLREW5.bat",
			"runLREW6.bat","runLREW7.bat","runLREW8.bat"))
# From the parameters given by DREAM, select out the parameter set for this chain.
this.par.set=pars[id,]
# Set the parameters in SWAT for this instance.
this.file<-model.in.file[id]
output.file<-output.rch.file[id]
write.table(this.par.set,this.file,
		quote=FALSE,row.names=model.in.rows,col.names=FALSE)
# Run the batch script that modifies the parameters and runs the model.
system(bfname[id])
# Read the model output.
Qsim.data<-read.table(output.file,skip=9)
Q1s<-Qsim.data[Qsim.data$V2==28,7]
Q2s<-Qsim.data[Qsim.data$V2==17,7] #The "17" here is the basin number of observed discharge data. THE
#				"6" is the column number in which the corresponding simulated discharge appears in
#				SWAT's"output.crh", and may very depending on options specified in SWAT's "file.cio".
# calculate the log likelihood. Here we are basing the likelihood function on the 
#  laplace distribution, with 0 as the location parameter and the maximum likelihood estimator for the scaler parameter
res1<-Q1o-Q1s
res2<-Q2o-Q2s
mean.res1<-mean(res1)
mean.abs.res1<-mean(abs(res1))
mean.res2<-mean(res2)
mean.abs.res2<-mean(abs(res2))
logp1<-sum(sapply(res1,function(res_k)
		log((1/(2*mean.abs.res1))*exp(-abs(res_k)/(mean.abs.res1)))))
logp2<-sum(sapply(res2,function(res_k)
		log((1/(2*mean.abs.res2))*exp(-abs(res_k)/(mean.abs.res2)))))
logp=logp1+logp2
		return(logp)
}
#
dream.object<-dream(
		FUN=callSWAT,
		pars=pars,
		func.type="logposterior.density",
		control=control)
summary(dream.object)
## Stop SNOW parallelisation
stopCluster(cl)
#
############################ END #########################################
