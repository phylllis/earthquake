# Extracting mainshock and largest aftershock

extract <- function(shocks){
  ms_data <- shocks[!is.na(shocks$ms),1:8]
  n.ms <- nrow(ms_data)
  as_data <- shocks[!is.na(shocks$as),]
  for (i in 1:n.ms){
    as.vector <- as_data[as_data$as==ms_data$ms[i],]$magnitude
    if(length(as.vector)==0) ms_data$as_max[i] <- NA
    if(length(as.vector)>0) ms_data$as_max[i] <- max(as.vector)
  }
  return(ms_data)
}