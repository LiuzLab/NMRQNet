
normalize=function(ppm,spec){
  ppm=rev(ppm)
  spec=rev(spec)
  area=0
  for(i in 1:(length(ppm)-1)){
    area=area+0.5*(spec[i]+spec[i+1])*(ppm[i+1]-ppm[i])
  }
  spec=spec/area
  return(spec)
}