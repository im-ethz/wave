dynamicIntegrateAUC <- function(survival.time, survival.status, start=NULL, marker, cutoffTime, event_code=NULL) {
   if(is.null(start))
      kmfit <- survfit(Surv(time=survival.time, event=survival.status) ~ 1)
   else
      kmfit <- survfit(Surv(time=start, time2=survival.time, event=survival.status) ~ 1)
   
   if(is.null(event_code)){
      mmm <- MeanRank( survival.time= survival.time, survival.status= survival.status, marker= marker, start=start )
   }
   else{
      mmm <- MeanRank_competing( survival.time= survival.time, survival.status= survival.status, marker= marker, start=start, event_code=event_code)
   }
   

   #Get overlap between survival function and mmm
   meanRanks <-  mmm$mean.rank[which(mmm$time <= cutoffTime)]
   survTimes <- mmm$time[mmm$time <= cutoffTime]
   timeMatch <- match(survTimes, kmfit$time)

   S_t <- kmfit$surv[timeMatch]

   #Calculate weights for c-index
   f_t <- c(1-S_t[1], (S_t[-length(S_t)] - S_t[-1]) )
   S_tao <- S_t[length(S_t)]
   weights <- (2*f_t*S_t)/sum(2*f_t*S_t)#(1-S_tao^2)
   
   return( sum(meanRanks * weights) ) #C-index
}

