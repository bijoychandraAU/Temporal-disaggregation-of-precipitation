# required libraries########  
library(dplyr, warn.conflicts = FALSE)
library(readxl)
library(lubridate)
library(readr)
library(data.table)
library(stringr)
library(MBC)

#1. name of observed file
obs_files="XX.csv"
obs_df <- lapply(obs_files, read_csv,col_types = cols())

#2. name of history file####
hist_files="XX_hist.csv"
hist_df <- lapply(hist_files, read_csv,col_names = FALSE,col_types = cols())
#3. name of future file########
fut_files="XX_future.csv"
rcp_df <- lapply(fut_files, read_csv,col_names = FALSE,col_types = cols())


# function for bias correcting the historical and future precipitation####
func<- function(obs1, hist1, rcp1,thres){
  rcp_ls=list()
  hist_ls=list()
  biased_hist_rcp <- list(hist_ls,rcp_ls)
  names(biased_hist_rcp) <- c("biased_hist", "biased_rcp")
  
  monthly_bias<- function(month){
    ############
    obs <- obs1
    colnames(obs)=c("time", "rain_in")
    
    obs=obs %>% 
      mutate(Precip_mm=25.4*rain_in) %>% 
      subset(select=-rain_in) %>% 
      #This is hourly rainfall data####
    #15 min to hourly observed rainfall
    mutate(datetime = floor_date(time, unit = "hour")) %>%
      group_by(datetime) %>%
      summarise(Rain_hour_mm = sum(Precip_mm, na.rm = TRUE))
    
    obs=obs %>% 
      mutate(mon=month(datetime)) %>% 
      filter(mon==month) %>% 
      subset(select=-(mon))
    
    
    #2. Read Historical data#########
    
    hist=hist1
    colnames(hist)=c("datetime","value")
    #2. Aggregation of hourly-historical precipitation to monthly (same steps as in 1)
    hist=hist %>%
      mutate(time=as.POSIXct(datetime,format="%Y-%m-%d%H:%M:%S",tz="UTC"),
             # Precip_his=value/24) %>% #for CANSEM and hadgem 
             Precip_his=value*3600) %>% # for GFDL and both MPI
      subset(select=-c(datetime,value))
    #round to nearUTC hour
    hist$time=format(round(hist$time, units="hours"), format="%Y-%m-%d %H:%M:%S")
    #remove duplicate values
    hist <- hist[!duplicated(hist[c('time')]),]
    #format change so  that there wont be NA in datetime
    hist$time=as.POSIXct(hist$time,format="%Y-%m-%d%H:%M:%S",tz="UTC")
    
    hist=hist %>% 
      mutate(mon=month(time)) %>% 
      filter(mon==month) %>% 
      subset(select=-(mon)) %>% 
      filter(Precip_his>0)
    #############
    threhold=sort(hist[["Precip_his"]])[thres*length(hist[["Precip_his"]])]
    
    
    hist[["Precip_his"]][hist[["Precip_his"]] <threhold] <- 0
    #Climate model which is for bias correction############
    
    hist=hist %>%  
      filter(Precip_his>0)
    
    
    rcp=rcp1
    colnames(rcp)=c("datetime","value")
    #3. Aggregation of hourly-future precipitation to monthly (same steps as in 1)#####
    rcp=rcp %>% 
      mutate(datetime=as.POSIXct(datetime,format="%Y-%m-%d%H:%M:%S",tz="UTC"),
             # Precip_fut=value/24) %>% #for CANSEM and hadgem 
             Precip_fut=value*3600) # for GFDL and both MPI
    
    rcp$datetime=format(round(rcp$datetime, units="hours"), format="%Y-%m-%d %H:%M:%S")
    #remove duplicate
    rcp <- rcp[!duplicated(rcp[c('datetime')]),]
    rcp$datetime=as.POSIXct(rcp$datetime,format="%Y-%m-%d%H:%M:%S",tz="UTC")
    
    rcp=rcp %>% 
      subset(select=-c(value)) %>% 
      mutate(mon=month(datetime)) %>% 
      filter(mon==month) %>% 
      subset(select=-(mon))%>% 
      filter(Precip_fut>0)

    ############
    rcp[["Precip_fut"]][rcp[["Precip_fut"]] <threhold] <- 0
    
    rcp=rcp %>%  
      filter(Precip_fut>0)
    
    ################
    
    obs_mat=as.matrix(obs$Rain_hour_mm)
    colnames(obs_mat) <- c("obs")
    hist_mat=as.matrix(hist$Precip_his)
    colnames(hist_mat) <- c("hist")
    rcp_mat=as.matrix(rcp$Precip_fut)
    colnames(rcp_mat) <- c("rcp")
    mylist <- list(obs_mat,hist_mat,rcp_mat,ratio.seq=TRUE)
    names(mylist) <- c("obs", "hist","rcp","ratio.seq")
    
    
    
    qdm.hist <- mylist$hist*0
    qdm.rcp <- mylist$rcp*0
    for(i in seq(ncol(mylist$hist))){
      
      fit.qdm <- QDM(o.c=mylist$obs[,i], m.c=mylist$hist[,i],
                     m.p=mylist$rcp[,i], ratio=mylist$ratio.seq[i],
                     # trace=mylist$trace1[i]
      )
      
      qdm.hist[,i] <- fit.qdm$mhat.c
      qdm.rcp[,i] <- fit.qdm$mhat.p
    }
    
    
    hist_bias=cbind(hist,qdm.hist)
    
    hist_bias=hist_bias %>% 
      subset(select=-(Precip_his))
    
    biased_hist_rcp$biased_hist[[i]]=hist_bias
    
    ##############
    rcp_bias=cbind(rcp,qdm.rcp)
    
    rcp_bias=rcp_bias %>% 
      subset(select=-(Precip_fut))
    
    biased_hist_rcp$biased_rcp[[i]]=rcp_bias
    ###########
    return(biased_hist_rcp)
  }

  month_ls=1:12
  batchSets <- split(month_ls, seq(length(month_ls)))
  bias_C=lapply(batchSets,monthly_bias)
  
  # unlist the list of list
  biased_unlist=do.call(c, unlist(bias_C, recursive=FALSE))
  
  # Historical part#######
  b_hist=biased_unlist[grep("biased_hist", names(biased_unlist))] 
  hist_correc=do.call("rbind", b_hist)
  hist_correc=hist_correc %>% 
    arrange(time)
  
  head(hist_correc)
  tail(hist_correc)
  # Future  part#######
  b_future=biased_unlist[grep("biased_rcp", names(biased_unlist))] 
  rcp_correc=do.call("rbind", b_future)
  rcp_correc=rcp_correc %>% 
    arrange(datetime)
  #set directory for storing the output
  write_csv(rcp_correc, file="biased_fut.csv")
  #set directory for storing the output
  write_csv(hist_correc, file="biased_hist.csv")
}
biased_out=parallel::mcmapply(function(obs1, hist1, rcp1,thres){
  return(func(obs1, hist1,rcp1,thres))
}, obs1=obs_df, hist1=hist_df,rcp1=rcp_df,thres=0.9,
mc.cores=1,SIMPLIFY = FALSE)





