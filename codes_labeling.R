label.shocks <- function(shocks){
  shocks <- cbind(shocks,'ms'=NA,'as'=NA,'fs'=NA)
  
  # Pre-specify the windows for defining aftershocks from Garden and Knopoff (1974)
  i.list <- 1:12 # grid number
  m.list <- seq(2.5,8,by=.5) # magnitude
  l.list <- c(19.5,22.5,26,30,35,40,47,54,61,70,81,94) # distance
  t.list <- c(6,11.5,22,42,83,155,290,510,790,915,960,985) # time elapsed
  
  # Pre-specify the window for defining foreshocks: 11 km, 7 days
  fs.l <- 11
  fs.t <- 7
  
  library(geosphere)
  # Sort the dataset.
  # First identify the largest shock in the dataset and label it as a mainshock.
  # Next identify all its foreshocks and aftershocks according to the above windows.
  # Then identify the next largest unlabelled shock in the data set and label it as a mainshock.  
  # Continue and repeat.
  
  shocks0 <- shocks[order(shocks$date,decreasing=F),]
  shocks0 <- shocks0[order(shocks0$magnitude,decreasing=T),]
  
  i <- 1
  k <- 1
  # while(nrow(shocks0)>0){
  while(shocks0$magnitude[1]>=5){
    ms.id <- shocks0$index[1]
    # LABEL MAIN SHOCK
    shocks[ms.id,]$ms <- i
    mainshock <- shocks[ms.id,] # mainshock entry
    # identify the aftershocks
    ms.mag <- mainshock$magnitude # mainshock magnitude
    ms.level <- floor((ms.mag-2)*2) # mainshock level in the window
    dist.ms <- function(ll){
      distm(ll,c(mainshock$longitude,mainshock$latitude))
    }
    as.l <- l.list[ms.level] # distance window for aftershocks
    as.t <- t.list[ms.level] # time window for aftershocks
    shocks.as <- shocks0[shocks0$date-mainshock$date>0 & shocks0$date-mainshock$date<=as.t*24*60*60,]
    dist.to.main.as <- apply(cbind(shocks.as$longitude,shocks.as$latitude),1,dist.ms)
    shocks.as2 <- shocks.as[dist.to.main.as <= as.l*1000,] # matrix of aftershocks
    # LABEL AFTER SHOCKS
    if (nrow(shocks.as2)>0) shocks[shocks.as2$index,]$as <- i
    # Identify the foreshocks
    shocks.fs <- shocks0[shocks0$date-mainshock$date<0 & shocks0$date-mainshock$date>=-fs.t*24*60*60,]
    dist.to.main.fs <- apply(cbind(shocks.fs$longitude,shocks.fs$latitude),1,dist.ms)
    shocks.fs2 <- shocks.fs[dist.to.main.fs <= fs.l*1000,] 
    # LABEL FORE SHOCKS
    if (nrow(shocks.fs2)>0) shocks[shocks.fs2$index,]$fs <- i
    # DELETE LABELED ENTRIES
    shocks0 <- shocks0[!(shocks0$index %in% c(ms.id,shocks.fs2$index,shocks.as2$index)),]
    print(i)
    i <- i+1
  }
  # re-index according to time
  new.index.list <- sum(!is.na(shocks$ms)):1
  old.index.list <- shocks$ms[!is.na(shocks$ms)]
  new.ms.index <- shocks$ms
  new.ms.index[!is.na(new.ms.index)] <- new.index.list
  index.change.fn <- function(x){
    if (is.na(x)==T){NA}
    else{sum(!is.na(shocks$ms))-which(old.index.list==x)+1}
  }
  new.as.index <- sapply(shocks$as,index.change.fn)
  new.fs.index <- sapply(shocks$fs,index.change.fn)
  shocks$ms <- new.ms.index
  shocks$as <- new.as.index
  shocks$fs <- new.fs.index
  return(shocks)
}