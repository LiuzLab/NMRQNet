

interpolation_fun=function(input_ppm,input_spectrum,length,min,max){
  l=length-length(input_spectrum)
  x=runif(l,min=min,max=max)
  update=approx(input_ppm,input_spectrum,xout=x)
  insert_ppm=update$x
  insert_intensity=update$y
  order=order(insert_ppm,decreasing=TRUE)
  insert_ppm=insert_ppm[order]
  insert_intensity=insert_intensity[order]
  update_data=as.data.frame(cbind(c(input_ppm,insert_ppm),c(input_spectrum,insert_intensity)))
  colnames(update_data)=c("ppm","intensity")
  update_data=update_data[order(update_data$ppm,decreasing=TRUE),]
  return(update_data$intensity)
}

