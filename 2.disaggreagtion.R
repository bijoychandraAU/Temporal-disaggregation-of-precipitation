# This R-code is for disaggregation of observed-15 min rainfall data.
# Here, the observed rainfall is recorded at 15-min interval.
# The data is aggregated back to hour. Then, it is disaggregated back to 15-min interval.

##Libraries####
suppressPackageStartupMessages({
  library(Hmisc) #for time lag
  library(tidyverse)
  library(dplyr, warn.conflicts = FALSE)
  library(readxl)
  library(lubridate)
  #MAE & RRSE used the following two libraries
  library(scorer)
  library(Metrics)
  library(plotly)
  library(beepr)
})

#Function for disaggregation of rainfall####

# name of observed csv files
file_name="XX.csv"
rain_15 <- function(j){
  df=read_csv(file_name,col_names = TRUE, col_types = cols())
  #rename and removed columns
  df=df %>% 
    mutate(Precip_mm=25.4*rain_inch) %>% 
    subset(select = -c(rain_inch))
  
  #convert inf to NAs
  is.na(df) <- sapply(df, is.infinite)
  #remove rows with NAs
  df =  na.omit(df)
  
  
  time1=df$datetime
  #it lagged time step by 1. example: 1,2,3 will become NA,1,2
  time2=Lag(df$datetime, shift = 1)
  # Creates eveent for all dataset having dry weather of 1-hr (it can be changed)
  df=df %>% 
    #dur=difference of time2 and time1 in hours
    mutate(dur=as.numeric(difftime(time1, time2, units = "hours"))) %>% 
    #count as an event only when duration>1 hours
    mutate(event=ifelse(dur>=1,1,0)) %>% 
    #add 1 to the previous event==1
    mutate(event1=event+c((event[-1] == 1) * 1, 0)) %>% 
    mutate(event=event1-event) %>% 
    subset(select = -event1)
  df =  na.omit(df)
  # add storm number and its length of each storm event
  df=df %>%
    #sum of each event till event=1
    group_by(Storm = as.integer(factor(lag(cumsum(event), default = 0)))) %>%
    #number of rows for each group
    mutate(strmlen = n()) 
  # Sub divided rainfall event of 1-hour i.e. when stormevent>1 hr are made into
  #1-hour interval
  #here 4 is used as it represents 1-hour
  MAX_GROUP = 4
  df=df %>% 
    group_by(Storm) %>% 
    #sub group for each group if nrow>4
    # fllor gives the numeric par of any number eg. floor(0.25)=0
    mutate(event_no = 
             paste(Storm, floor((row_number()-1)/ (MAX_GROUP)), sep = "-")) %>% 
    ungroup()
  # Calculate depth of each event after sub-divided
  df=df %>% 
    group_by(event_no) %>% 
    #depth of rainfall events with 1-hour
    mutate(Depth = sum(Precip_mm))
  df=df %>% 
    subset(select=-c(event,Storm,strmlen))%>%
    mutate(month = month(datetime))
  df=df[order(df$Depth),]
  # End of event database for 1-hr####
  #2.event for each month in ascending order####
  df_month=function(mon){
    month_eventlist = list()
    df_mon=subset(df, month==mon)
    df_mon=df_mon[order(df_mon$Depth),]
    month_eventlist[[j]] <- df_mon
    return(month_eventlist)
  }
  
  list_final=mapply (function(mon) df_month(mon),mon=01:12, SIMPLIFY = FALSE)
  list_final1 <- unlist(list_final, recursive = FALSE)

  event_station<- do.call("rbind", list_final1)
  epsilon=min(event_station$Depth)
  
  #3. Aggregation of observed -15min to hourly data####
  # read observed discontinuous 15-min rainfall
  df_hour1<- subset(df, select = c(datetime,Precip_mm ))
  df_hour1=df_hour1[order(df_hour1$datetime),]
  # aggregation of 15-min continuous data to hourly data####
  df_hour1= df_hour1 %>%
    mutate(datetime = floor_date(datetime, unit = "hour")) %>%
    group_by(datetime) %>%
    summarise(Rain_hour_mm = sum(Precip_mm, na.rm = TRUE)) %>% 
    mutate(month = month(datetime))
  df_hour=df_hour1[df_hour1$Rain_hour_mm != 0, ]

  # 4.Disaggregation of hourly to quarter-hourly precipitation#####
  rain15 <- function(x,y,mon){
    rain15_station=list()
    #create CDF
    df_mont=df %>% 
      filter(month==mon)
    
    cdf1=df_mont[,c("Depth")]
    cdf2 <- ecdf(cdf1$Depth)
    
    epsilon=min(cdf1$Depth)
    
    Ddt=x
    if (Ddt>epsilon) {
      mat=list()
      y=y
      Dt=Ddt
      i<-0
      while (TRUE){
        i<-i+1
        a=cdf2(Dt)
        u1 <- runif(1, 0, a)
        d1=quantile(cdf2, u1)
        ###filtering with raindepth=d1
        event_rain=df_mont %>% 
          # find the value that is closest to d1
          subset(Depth == max(Depth[Depth<= d1])) 
        #random selection of rain event==d1
        event_rand=event_rain[sample(nrow(event_rain), 1),]
        D1=event_rain%>% 
          filter(event_no==event_rand$event_no) %>% 
          subset(select=c("datetime","Precip_mm"))
        #########
        y=as.data.frame(y[1])
        colnames(y)=c("datetime")
        df_15=y%>%
          add_row(datetime = as_datetime(y$datetime)+minutes(15)) %>%
          add_row(datetime = as_datetime(y$datetime)+minutes(30)) %>%
          add_row(datetime = as_datetime(y$datetime)+minutes(45))%>% 
          mutate(Precip_mm=0)
        #Random insert of data and add the final value####
        random_row1 <- sample.int(n=nrow(df_15),
                                  size=nrow(D1),
                                  replace = FALSE)
        df_15$Precip_mm[random_row1]<-df_15$Precip_mm[random_row1] + D1$Precip_mm
        mat[[i]]<-df_15
        Dt<-Dt-sum(df_15[,"Precip_mm"])
        if (Dt <= epsilon){break}
      }
      ###sum all rain for each group of list
      rain_15f= mat %>%
        bind_rows() %>%
        group_by(datetime ) %>%
        summarise_all(sum)
      
      d_remain=Ddt-sum(rain_15f[,2])
      
      rand_row <- sample.int(n=nrow(rain_15f),
                             size=length(d_remain),
                             replace = FALSE)
      rain_15f$Precip_mm[rand_row]<-rain_15f$Precip_mm[rand_row] + d_remain
      
    } else {
      #random insertion of remaining rain < minimum rain
      y=as.data.frame(y[1])
      colnames(y)=c("datetime")
      rain_15f=y%>%
        add_row(datetime = as_datetime(y$datetime)+minutes(15)) %>%
        add_row(datetime = as_datetime(y$datetime)+minutes(30)) %>%
        add_row(datetime = as_datetime(y$datetime)+minutes(45))%>% 
        mutate(Precip_mm=0)
      
      rand_row <- sample.int(n=nrow(rain_15f),
                             size=1,
                             replace = FALSE)
      rain_15f$Precip_mm[rand_row]<-rain_15f$Precip_mm[rand_row] + Ddt
    }
    return(rain_15f) 
  }
  
  rain_15times=mapply (function(x,y,mon) rain15(x,y,mon),
                       x=df_hour$Rain_hour_mm,y=df_hour$datetime, 
                       mon=df_hour$month,SIMPLIFY = FALSE)
  
  rain_times <- do.call("rbind", rain_15times)
  
  start=df_hour1$datetime[1]
  end=df_hour1$datetime[nrow(df_hour1)]
  DateTime=(seq(start,end+60*60*1, by="15 min"))
  
  simu_15min=as_tibble(DateTime)
  colnames(simu_15min)="datetime"
  simu_15min=simu_15min[-nrow(simu_15min),]
  #check for dim as discontinuous df_hour
  #this is the continuous simulated/disaggregated 15-min rainfall
  simu_15min=left_join(simu_15min, rain_times, by = 'datetime') %>%
    mutate(Precip_mm = replace_na(Precip_mm, 0)) %>% 
    mutate(Precip_mm = round(Precip_mm, 2))
  write_csv(simu_15min, file="obs_15.csv")
  print(j)
}
output15<- replicate(1,t(mapply(function(j) rain_15(j),
                                 j=file_name)))

