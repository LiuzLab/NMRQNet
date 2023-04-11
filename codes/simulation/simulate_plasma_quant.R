##################Simulation pipeline for plasma samples############################

library(parallel)
source("./main_functions/integral.R")
source("./main_functions/generate_lorent.R")

##Initialize the chemical shift values
ppm_region=seq(0.7,4.1,length.out = 12000)
ppm_step=(4.1-0.7)/12000

###Load 38 metabolites peak signatures
metabolite_lib_info=read.csv("reference_library_38.csv",sep=",")
metabolite_lib_info=metabolite_lib_info[,-1]

#Load lipoprotein clusters
cus_lipo_shift=read.csv("cus_lipo_data_shift_intergral.csv")
cus_lipo_shift=cus_lipo_shift[,-1]
cus_lipo_shift=as.matrix(cus_lipo_shift)
colnames(cus_lipo_shift)=NULL

##Load DSS data
hp_DSS_info=read.csv("real_DSS_info_integral.csv",sep=",",header = TRUE)
hp_DSS_info=hp_DSS_info[,-1]

##Load anticoagulants spectra
#EDTA as anticoagulant
EDTA_dss=read.csv("real_EDTA_dss_info_integral.csv",sep=",")
EDTA_dss=EDTA_dss[,-1]
EDTA_anti=EDTA_dss$intensity

##Sodium_heparin as anticoagulant
sheparin_dss=read.csv("real_sHeparin_dss_info_integral.csv",sep=",")
sheparin_dss=sheparin_dss[,-1]
sheparin_anti=sheparin_dss$intensity

## concentration table
concentration_table=rbind(c(0.73,0.52),c(0.55,0.01),c(0.82,29.94),
                          c(1.08,0.20),c(0.50,5.00),c(1.69,10.67),c(0.52,0.20))
concentration_table=as.data.frame(concentration_table)
rownames(concentration_table)=c("organoids","human_plasma","cerebellum","cells",
                                "amniotic_fluid","fly_brain","cell_media")
colnames(concentration_table)=c("shape","rate")

d=2 #use the fitted curve from human plasma annotated results
shape=concentration_table$shape[d]
rate=concentration_table$rate[d]

dir="./quant_concen_train/" #saved directory

sim_fun=function(n,shape,rate,ppm_region,ppm_step,metabolite_lib_info,
                 EDTA_anti,sheparin_anti,hp_DSS_info){
  
  num=sample(10:38,1) #mix at least 10 metabolites in each synthesized spectrum
  shift=runif(1,-0.03,0.03) #positional variations
  step=floor(shift/ppm_step)
  
  compound=sample(1:38,num) #metabolite candidates in the mixture
  meta_quant=rep(0,38) #Initialize concentrations for metabolites and other components
  lipo_quant=rep(0,24)
  anti_quant=rep(0,3)
  
  ##generate whole spectrum for 38 metabolites
  tem_library=NULL
  for(k in 1:length(compound)){
    
    #compound-level key, use key to regularize the pattern dynamic variations
    key=sample(1:3,1)
    
    index=compound[k]
    meta=as.character(unique(metabolite_lib_info$Name)[index])
    list=which(metabolite_lib_info$Name==meta) #extract the peak signatures from the library based on the metabolite name
    
    spec=rep(0,12000)
    for(j in 1:length(list)){
      if((metabolite_lib_info$loc[list[j]]<4.1) & 
         (metabolite_lib_info$loc[list[j]]>0.7) & (key==1)){
        scaler=450 #simulate sharp peaks (high resolution)
        spec=spec+generate_lorent(ppm_region,metabolite_lib_info[list[j],],scaler=scaler) }
      else if((metabolite_lib_info$loc[list[j]]<4.1) & (metabolite_lib_info$loc[list[j]]>0.7) & (key==2)){
        scaler=runif(1,min=200,max=300) #simulate wide peaks or with overlappings
        spec=spec+generate_lorent(ppm_region,metabolite_lib_info[list[j],],scaler=scaler)  }
      else if((metabolite_lib_info$loc[list[j]]<4.1) & (metabolite_lib_info$loc[list[j]]>0.7) & (key==3)){
        scaler=runif(1,min=200,max=1000) #simulate dynamic peak scales
        spec=spec+generate_lorent(ppm_region,metabolite_lib_info[list[j],],scaler=scaler) }
    }
    tem_library=rbind(tem_library,spec)
  }
  ##for all generated metabolite spectra, scale to have the area under the curve equal to 1
  metabolite_spec=t(apply(tem_library,1,function(x){return(x/integral(ppm_region,x))}))
  
  ##randomly select the anticoagulant in each simulated spectrum
  anti_list=c(1,2,3,3)
  anti=sample(anti_list,1)
  if(anti==1){
    anti_compound=EDTA_anti
  }else if(anti==2){
    anti_compound=sheparin_anti
  }else{
    anti_compound=hp_DSS_info$intensity
  }
  
  concen_num=num+25 #metabolites + 1 anticoagulant + 24 lipoprotein clusters
  concentration=rgamma(concen_num,shape=shape,rate=rate)
  
  anti_concen=sample(1:concen_num,1)
  lipo_concen=sample((1:concen_num)[-anti_concen],24)
  compound_concen=(1:concen_num)[-c(anti_concen,lipo_concen)]
  
  input=colSums(c(concentration[compound_concen])*metabolite_spec) 
  input=input+colSums(c(concentration[lipo_concen])*cus_lipo_shift)
  input=input+concentration[anti_concen]*anti_compound #combine different components and their concentrations
  
  #embed positional variations
  if(step<0){
    input[(1:(12000-abs(step)))]=input[(abs(step)+1):12000]
    input[((12000-abs(step)+1):12000)]=rep(0,abs(step))
  }
  if(step>0){
    input[(step+1):12000]=input[(1:(12000-step))]
    input[(1:step)]=rep(0,step)
  }
  
  meta_quant[compound]= concentration[compound_concen]
  lipo_quant=concentration[lipo_concen]
  if(anti==1){
    anti_quant[3]=concentration[anti_concen]
  }else if(anti==2){
    anti_quant[2]=concentration[anti_concen]
  }else{
    anti_quant[1]=concentration[anti_concen]
  }
  
  quant_output=c(meta_quant,lipo_quant,anti_quant)
  
  quant_output=quant_output/integral(ppm_region,input)
  input=input/integral(ppm_region,input)
  
  #save out the simulated mixture spectrum and the concentrations for different components
  input_ID=paste("raw_input",n,sep="_")
  input_ID=paste(input_ID,"csv",sep=".")
  input_path=paste(dir,input_ID,sep="")
  
  quant_concen_ID=paste("quant_concen",n,sep="_")
  quant_concen_ID=paste(quant_concen_ID,"csv",sep=".")
  quant_concen_path=paste(dir,quant_concen_ID,sep="")
  
  write.csv(input,input_path)
  write.csv(quant_output,quant_concen_path)
  
}

mc <- getOption("mc.cores", 20)

mclapply(1:20000, FUN=sim_fun,shape,rate,ppm_region,ppm_step,
         metabolite_lib_info,
         EDTA_anti,sheparin_anti,hp_DSS_info, mc.cores = mc)

