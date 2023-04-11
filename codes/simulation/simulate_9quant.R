####### Simulate the training set for nine-metabolite quantification study####
library(parallel)
source("./main_functions/integral.R")

#setup ppm
ppm=seq(-1,11,length.out = 40000)
ppm_step=12/40000

#9 metabolites normalized synthesized spectra
syn_9_norm=read.csv("chenomx_9meta_alignreal.csv",sep=",")
rownames(syn_9_norm)=as.character(syn_9_norm[,1])
syn_9_norm=syn_9_norm[,-1]
colnames(syn_9_norm)=NULL

#9 metabolites normalized real spectra (with DSS)
real_9_norm=read.csv("real_9meta_norm_align.csv",sep=",")
rownames(real_9_norm)=as.character(real_9_norm[,1])
real_9_norm=real_9_norm[,-1]
colnames(real_9_norm)=NULL

##Estimated dss/meta area ratio in real 9-metabolite NMR spectra
dss_meta_ratio=read.csv("dss_meta_ratio_renew.csv",sep=",")
colnames(dss_meta_ratio)=c("meta_name","DSS_ratio","meta_ratio")

##Fitted gamma distributions from annotated concentrations
concentration_table=rbind(c(0.73,0.52),
                          c(0.55,0.01),
                          c(0.82,29.94),
                          c(1.08,0.20),
                          c(0.50,5.00),
                          c(1.69,10.67),
                          c(0.52,0.20))
concentration_table=as.data.frame(concentration_table)
rownames(concentration_table)=c("organoids","human_plasma","cerebellum","cells","amniotic_fluid","fly_brain","cell_media")
colnames(concentration_table)=c("shape","rate")


dir="./train_9quant/" #saved path
d=sample(1:7,1) #random choose one gamma function
shape=concentration_table$shape[d]
rate=concentration_table$rate[d]

sim_fun=function(n,shape,rate,ppm,ppm_step,syn_9_norm,real_9_norm,
                 dss_meta_ratio){

  num=sample(1:9,1) #number of metabolites
  shift=runif(1,-0.03,0.03) #bias in chemical shift
  step=floor(shift/ppm_step) # resulted shift steps
  
  compound=sample(1:9,num) ##metabolite candidates within the mixture
  concen=rep(0,10) ##initialize concentrations for ten components
  
  real_index=sample(1:length(compound),1) #random choose one metabolite and mix its
  #real NMR spectrum within the mixture
  real_compound=compound[real_index]
  
  concentration=rgamma(num,shape=shape,rate=rate) #generate concentrations from fitted gamma distr.
  concen[real_compound]=concentration[1]
  
  if(num>1){
    syn_compound=compound[-real_index] #for the rest of metabolites, their synthesized spectra are mixed
    syn_spec=syn_9_norm[syn_compound,]
    
    ##leave the first concentration for the real metabolite
    ## since DSS does not shift, so we only shift synthesized spectra
    syn_spectra=colSums(c(concentration[-1])*as.matrix(syn_spec))
    concen[syn_compound]=c(concentration[-1])
    
    if(step<0){
      syn_spectra[(1:(40000-abs(step)))]=syn_spectra[(abs(step)+1):40000]
      syn_spectra[((40000-abs(step)+1):40000)]=rep(0,abs(step))
    }
    if(step>0){
      syn_spectra[(step+1):40000]=syn_spectra[(1:(40000-step))]
      syn_spectra[(1:step)]=rep(0,step)
    }
    
    total_spectra=colSums(rbind(concentration[1]*as.numeric(real_9_norm[real_compound,]),
                                syn_spectra))
    
  }else{
    total_spectra=concentration[1]*as.numeric(real_9_norm[real_compound,])
  }
  
  ##estimate the DSS and metabolite areas within the real spectrum
  ratio_DSS=dss_meta_ratio$DSS_ratio[real_compound]
  ratio_meta=dss_meta_ratio$meta_ratio[real_compound]
  
  #10 is the index for DSS in the output
  concen[10]=concen[real_compound]*ratio_DSS
  concen[real_compound]=concen[real_compound]*ratio_meta
 
  concen=concen/integral(ppm,total_spectra)
  total_spectra=total_spectra/integral(ppm,total_spectra)
  total_spectra[which(total_spectra<0)]=0
  
  set_ID=paste("training_set",n,sep="_")
  concen_ID=paste("training_concen",n,sep="_")
  set_ID=paste(set_ID,"csv",sep=".")
  concen_ID=paste(concen_ID,"csv",sep=".")
  set_path=paste(dir,set_ID,sep="")
  concen_path=paste(dir,concen_ID,sep="")
  
  write.csv(total_spectra,set_path)
  write.csv(concen,concen_path)
}

mc <- getOption("mc.cores", 20)

mclapply(1:10000, FUN=sim_fun,
         shape=shape,rate=rate,ppm=ppm,
         ppm_step=ppm_step,
         syn_9_norm=syn_9_norm,real_9_norm=real_9_norm,
         dss_meta_ratio=dss_meta_ratio, mc.cores = mc)

