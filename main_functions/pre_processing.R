pre_processing_fun=function(ppm,intensity,a,b){
  mat=intensity[which(ppm>=a&ppm<=b)]
  mat=as.matrix(mat)
  mat=t(mat)
  rownames(mat)=c("Ref")
  colnames(mat)=as.character(ppm[which(ppm>=a&ppm<=b)])
  spectra_cor=BaselineCorrection(mat)
  spectra_cor[1,which(spectra_cor[1,]<0)]=0
  return(spectra_cor[1,])
  
 # plot_data_frame=cbind(ppm[which(ppm>=a&ppm<=b)],spectra_cor[1,])
#  colnames(plot_data_frame)=c("ppm","intensity")
#  plot_data_frame=as.data.frame(plot_data_frame)
#  return(plot_data_frame)
  
}