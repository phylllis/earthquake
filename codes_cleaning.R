clean.shocks <- function(shocks){
  # cleaning the time variable
  shocks <- shocks[,c(2,2,2,2,2,5,6,7,11,12)]
  colnames(shocks) <- c('date','year','month','day','time','latitude','longitude','depth','type','magnitude')
  shocks$date <- as.character(shocks$date)
  timepoint <- matrix(unlist(strsplit(shocks$date,' ')),ncol=2,byrow=T)
  timepoint_ymd <- matrix(unlist(strsplit(timepoint[,1],'/')),ncol=3,byrow=T)
  shocks$time <- timepoint[,2]
  shocks$year <- timepoint_ymd[,3]
  shocks$month <- timepoint_ymd[,1]
  shocks$day <- timepoint_ymd[,2]
  remove(timepoint_ymd,timepoint)
  # The time points only recorded the last two digits of the year.  The first two digits will be manually added.  
  century.pt <- 1319
  century.prefix <- c(rep(20,century.pt),rep(19,nrow(shocks)-century.pt))
  shocks$year <- paste(century.prefix,shocks$year,sep='')
  shocks$date <- strptime(paste(shocks$year,'-',shocks$month,'-',shocks$day,' ',shocks$time,sep=''),'%Y-%m-%d %H:%M')
  shocks <- shocks[,-(2:5)]
  remove(century.pt,century.prefix)
  # cleaning the other variables
  shocks$latitude <- as.numeric(as.character(shocks$latitude))/10
  shocks$longitude <- as.numeric(as.character(shocks$longitude))/10
  shocks$depth <- as.numeric(as.character(shocks$depth))
  shocks$magnitude <- as.numeric(as.character(shocks$magnitude))
  index <- 1:nrow(shocks)
  shocks <- cbind(index,shocks)
  remove(index)
  return(shocks)
}