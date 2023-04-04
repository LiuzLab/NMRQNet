generate_lorent=function(ppm,input,scaler){
  A=input$height
  W=input$width/scaler
  loc=input$loc
  
  output=A/(1+((ppm-loc)/(W/2))^2)
  return(output)
}