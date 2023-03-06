## max_pooling return both down-resolution spectrum and the chosen index

max_pool=function(spectrum,window_size){
  down_spec=NULL
  ppm_down=NULL
  step=length(spectrum)/window_size
  for(i in c(1:step)){
    down_spec=c(down_spec,max(spectrum[c(((i-1)*window_size+1):(i*window_size))]))
    ppm_down=c(ppm_down,window_size*(i-1)+which.max(spectrum[c(((i-1)*window_size+1):(i*window_size))]))
  }
  return(cbind(down_spec,ppm_down))
}