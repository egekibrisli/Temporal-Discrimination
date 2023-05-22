# Get time productions

library(RcppCNPy)
library(R.utils)
library(scales)
library(LambertW)
library(ggplot2)
library(ggExtra)
library(mgcv)
library(lme4)
library(lmerTest)
library(dplyr)
library(itsadug)
library(plyr)
library(stringr)


#=======================================================================================
# UTILS & functions
#=======================================================================================
make_trial_index <- function(d) {
	trial_start_stamp <- match(d$V5, 'StartPoke') + match(d$V3, 'Entry') - 1
	sum(trial_start_stamp, na.rm=T)
	time_stamps <- which(trial_start_stamp==1)
	time_stamps <- c(1,time_stamps, length(trial_start_stamp))
	trial_index <- c()
	for (i in 2:length(time_stamps)) {
		trial_index <- c(trial_index, rep(i-1, time_stamps[i]-(time_stamps[i-1])))
	}
	trial_index <- c(trial_index, i)
	d <- cbind(d, trial_index)
	d <- subset(d, trial_index!=1)
	return(d)
}


#=======================================================================================
# UTILS & functions
#=======================================================================================
dat <- data.frame()
for (ss in c(1:15)) { # sessions number
  dir <- paste('C:/Users/PC/Dropbox/eloras_rats/data/step6/S', toString(ss),'/', sep='')
	ls <- list.files(dir)
	csv_list <- c()
	for (i in c(1:length(ls))) {
		if (grepl('.csv', ls[i])) {
			csv_list <- c(csv_list, ls[i])
		}}
	for (n in c(1:length(csv_list))) {
		d <- read.table(paste(dir, csv_list[n], sep=''), header=F, nrows=4, sep=',', fill=T)
		rat_nr <- toString(d$V4[1])
		group_width <- substr(toString(d$V4[2]), start=1, stop=4)
		group_dir <- substr(toString(d$V4[2]), start=6, stop=7)
		print(rat_nr)
		print(group_width)
		print(group_dir)
		d <- read.table(paste(dir, csv_list[n], sep=''), header=F, skip=7, sep=',', fill=T)
		d <- d[-c((nrow(d)-3):nrow(d)), ] # remove last 4 line which is trash
		d$time_stamp <- d$V2 #+ d$V3/10^3
		d <- make_trial_index(d)
		d_len <- dim(d)[1]
		d <- d[-c(nrow(d)), ]
		d <- cbind(d, rep(ss,length(d_len)), rep(group_width, length(d_len)), rep(group_dir, length(d_len)), rep(rat_nr,length(d_len)))
		dat <- rbind(dat, d)
	}}

colnames(dat) <- c('V1', 'V2', 'V3', 'V4', 'V5','V6','V7','time_stamp', 'trial_index', 'session', 'group_width', 'group_dir', 'rat_nr')

dat$V6 <- as.numeric(as.character(dat$V6))

#=======================================================================================
# UTILS & functions
#=======================================================================================


all_data <- ddply(dat, .(rat_nr, session, trial_index, group_width, group_dir), function(z){ #int_offset
	#z <- subset(dat, rat_nr==9 & session==1 & trial_index==10) z is for debugging
	tmp1 <- subset(z, V5=='StartPoke' & V3=='Entry')
	tmp2 <- subset(z, V5=='StartPort' & V3=='Input')

	start_light = tmp1$time_stamp[1] # nose poke start
	init_entry = tmp2$time_stamp[1] # initial entry in start port

	cond=0
	if (dim(subset(z, V5=='Duration_1s' & V3=='Exit'))[1]>0) {cond=1 # find out condition
	} else if (dim(subset(z, V5=='Duration_2s' & V3=='Exit'))[1]>0) {cond=2
	} else if (dim(subset(z, V5=='Duration_4s' & V3=='Exit'))[1]>0) {cond=4}


	int_offset = subset(z, grepl('Duration_', V5) & V3=='Exit')$time_stamp[1]
	first_check = subset(z, grepl('NosePoke', V5) & time_stamp>int_offset)$V4[1] # check first state

	first_check_lat = subset(z, grepl('NosePoke', V5) & time_stamp>int_offset)$time_stamp[1]
	
	duration_stamps <- list()
	temp_durations <- list()
	for (i in 1:nrow(dat)) {
	  # check if the V5 value contains "Duration_" and V3 value is "Exit"
	  if (grepl("Duration_", dat[i, "V5"]) & dat[i, "V3"] == "Exit") {
	    # if so, add the timestamp to the temp_durations list
	    temp_durations[[length(temp_durations) + 1]] <- dat[i, "V2"]
	    
	    # check the next 7 rows for another "Duration_" value
	    for (j in i+1:min(i+7, nrow(dat))) {
	      if (grepl("Duration_", dat[j, "V5"]) & dat[i, "V3"] == "Exit") {
	        # if another "Duration_" value is found, clear the temp_durations list
	        temp_durations <- list()
	        # add the new timestamp to the temp_durations list
	        temp_durations[[length(temp_durations) + 1]] <- dat[j, "V2"]
	        # move to the next iteration of the outer loop
	        break
	      } else if (j == min(i+7, nrow(dat))) {
	        # if no other "Duration_" value is found within the next 7 rows, add the temp_durations
	        # list to the duration_stamps list and clear the temp_durations list
	        duration_stamps[[length(duration_stamps) + 1]] <- temp_durations
	        temp_durations <- list()
	      }
	    }
	  }
	}
	response_nose = subset(dat, grepl('NosePokeChoice_', V5) & V3 == 'Exit')$time_stamp

	
	true_duration_stamps <- lapply(duration_stamps, function(x) {
	  if(length(x) == 2) {
	    x <- x[-1]
	  }
	  return(x)
	})
	unlisted <- lapply(true_duration_stamps, function(x) {
	  if(is.list(x) && length(x) == 1) {
	    x <- unlist(x)
	  }
	  return(x)
	})
	
	response_latency <- mapply(function(x, y) {
	  x - y
	}, response_nose, unlisted)
	

	return(c(cond, start_light, init_entry, int_offset,first_check, first_check_lat, response_latency))
	}, .progress='text')
colnames(all_data)[c(6:12)] <- c('cond','start_light', 'init_entry', 'int_offset', 'first_check', 'first_check_lat','response_latency')


# wide group
all_data$correct_entry <-
	ifelse(all_data$group_width=='WIDE' & all_data$group_dir=='LM' & all_data$cond==1 & all_data$first_check==3, 1,
	ifelse(all_data$group_width=='WIDE' & all_data$group_dir=='LM' & all_data$cond==2 & all_data$first_check==5, 1,
	ifelse(all_data$group_width=='WIDE' & all_data$group_dir=='LM' & all_data$cond==4 & all_data$first_check==1, 1,
	ifelse(all_data$group_width=='WIDE' & all_data$group_dir=='ML' & all_data$cond==1 & all_data$first_check==3, 1,
	ifelse(all_data$group_width=='WIDE' & all_data$group_dir=='ML' & all_data$cond==2 & all_data$first_check==1, 1,
	ifelse(all_data$group_width=='WIDE' & all_data$group_dir=='ML' & all_data$cond==4 & all_data$first_check==5, 1,
	ifelse(all_data$group_width=='NARO' & all_data$group_dir=='LM' & all_data$cond==1 & all_data$first_check==3, 1,
	ifelse(all_data$group_width=='NARO' & all_data$group_dir=='LM' & all_data$cond==2 & all_data$first_check==4, 1,
	ifelse(all_data$group_width=='NARO' & all_data$group_dir=='LM' & all_data$cond==4 & all_data$first_check==2, 1,
	ifelse(all_data$group_width=='NARO' & all_data$group_dir=='ML' & all_data$cond==1 & all_data$first_check==3, 1,
	ifelse(all_data$group_width=='NARO' & all_data$group_dir=='ML' & all_data$cond==2 & all_data$first_check==2, 1,
	ifelse(all_data$group_width=='NARO' & all_data$group_dir=='ML' & all_data$cond==4 & all_data$first_check==4, 1, 0))))))))))))


write.table(all_data, 'C:/Users/PC/Dropbox/eloras_rats/data/all_data_step6.csv')













# #=======================================================================================
# # Grab all of the  data in one data frame first.
# #=======================================================================================
# dat <- data.frame()
# for (b in c(1:40)) {
# 	dir <- paste('/Users/tadeusz/Dropbox/Metacognition rat orsay DATA/metacognition_rat_orsay/data_testing/data_testing_step5_ss', toString(b),'/', sep='')
# 	ls <- list.files(dir)
# 	csv_list <- c()
# 	for (i in c(1:length(ls))) {
# 		if (grepl('.csv', ls[i])) {
# 			csv_list <- c(csv_list, ls[i])
# 		}
# 	}
# 	for (n in c(1:length(csv_list))) {
# 		d <- read.table(paste(dir, csv_list[n], sep=''), header=F, nrows=2, sep=',', fill=T)
# 		rat_nr <- toString(d$V4[1])
# 		group <- substr(toString(d$V4[2]), start=4, stop=7)
# 		if (substr(toString(d$V4[2]), start=4, stop=4) == 'P') {
# 			interval <- as.numeric(substr(toString(d$V4[2]), start=10, stop=11)) / 10
# 			} else {
# 			interval <- as.numeric(substr(toString(d$V4[2]), start=9, stop=10)) / 10
# 		}
# 		print(group)
# 		print(rat_nr)
# 		d <- read.table(paste(dir, csv_list[n], sep=''), header=F, skip=7, sep=',', fill=T)
# 		d <- d[-c((nrow(d)-3):nrow(d)), ] # remove last 4 line which is trash
# 		d$time_stamp <- d$V2 + d$V3/10^3
# 		d <- make_trial_index(d)
# 		d_len <- dim(d)[1]
# 		d <- d[-c(nrow(d)), ] # last line TO DO should be fixed in make_trial_index.R
# 		d <- cbind(d, rep(b,length(d_len)), rep(group, length(d_len)), rep(rat_nr,length(d_len)))
# 		dat <- rbind(dat, d)
# 	}
# }
# colnames(dat) <- c('V1', 'V2', 'V3', 'V4', 'V5','V6','V7','V8','V9','time_stamp', 'trial_index', 'session', 'group', 	'rat_nr')
# dat$rat_nr <- as.numeric(as.character(dat$rat_nr))
# # TO DO: Remove col 'V2, V3, V9, V8'
#
#
# #=======================================================================================
# # Debug >>>
# #z <- subset(dat, rat_nr==1 & trial_index==4 & session==6)
#
#
# all_data <- ddply(dat, .(rat_nr, group, session, trial_index), function(z){
# 	rts<-0
# 	pellet<-0
# 	pellet_test<-0
# 	if (z$group[1]=='HOLD') {
# 		tmp1 <- subset(z, V6=='Lever' & V4=='Input')
# 		tmp2 <- subset(z, V6=='Lever off' & V4=='Input')
# 		tp <- tmp2$time_stamp - tmp1$time_stamp
# 	} else if (z$group[1]=='PRES') {
# 		tmp <- subset(z, V6=='Lever' & V4=='Input')
# 		if (dim(tmp)[1]!=2) {tp<-0 # It catches last trial when program stoped due to time limit
# 			} else {
# 			tmp1 <- tmp[seq(1, dim(tmp)[1], 2), ]
# 			tmp2 <- tmp[seq(2, dim(tmp)[1], 2), ]
# 			tp <- tmp2$time_stamp - tmp1$time_stamp
# 		}
# 	}
# 	if (length(tp)>1) {tp<-0} else if (length(tp)==0) {tp<-0}
# 	# EXTRACT port side
# 	if (sum(z$V6=='FDR 1 0p Light' & z$V4=='Entry')) {reward_port=1
# 		} else if (sum(z$V6=='FDR 1 1p Light' & z$V4=='Entry')) {reward_port=1
# 		} else if (sum(z$V6=='FDR 1 2p Light' & z$V4=='Entry')) {reward_port=1
# 		} else if (sum(z$V6=='FDR 2 0p  Light' & z$V4=='Entry')) {reward_port=2
# 		} else if (sum(z$V6=='FDR 2 1p Light' & z$V4=='Entry')) {reward_port=2
# 		} else if (sum(z$V6=='FDR 2 2p Light' & z$V4=='Entry')) {reward_port=2
# 		} else if (sum(z$V6=='FDR 1&2 2p Light' & z$V4=='Entry')) {reward_port=32 # ?????
# 		} else if (sum(z$V6=='FDR 1&2 1p Light' & z$V4=='Entry')) {reward_port=31
# 		} else {reward_port=0}
# 	# EXTRACT RT: left and right, if no rt than assign 0.
# 	subset_to_eval <- subset(z, grepl('Light',V6) & grepl('Entry',V4))
# 	if (!empty(subset_to_eval)) {
# 		rts <- subset(z, grepl('Port',V6))$time_stamp - subset_to_eval$time_stamp
# 		# EXTRACT pellet number on a given trial:
# 		switch(as.numeric(substr(as.character(subset_to_eval$V6), start=7, stop=7))+1,
# 		{pellet=0}, {pellet=1}, {pellet=2})
# 		# EXTRACT pellet number on a given trial in "1&2" test trials:
# 		if (!empty(subset(z, grepl('FDR 1&2 1p Pellet',V6) & grepl('Entry',V4)))) {pellet_test=1
# 			} else if (!empty(subset(z, grepl('FDR 1&2 2p Pellet',V6) & grepl('Entry',V4)))) {pellet_test=2
# 			} else if (empty(subset(z, grepl('FDR 1&2 2p Pellet',V6) & grepl('Entry',V4))) & empty(subset(z, grepl('FDR 1&2 2p Pellet',V6) & grepl('Entry',V4)))) {pellet_test=0}
# 	}
# 	# Rewritten to 3 lines + pellet extraction -->
# 	# if (sum(z$V6=='FDR 1 0p Light' & z$V4=='Entry')) {
# 	# 	rts <- subset(z, grepl('Port',V6))$time_stamp - subset(z, V6=='FDR 1 0p Light' & V4=='Entry')$time_stamp
# 	# }	else if (sum(z$V6=='FDR 1 1p Light' & z$V4=='Entry')) {
# 	# 	rts <-	subset(z, grepl('Port',V6))$time_stamp - subset(z, V6=='FDR 1 1p Light' & V4=='Entry')$time_stamp
# 	# }	else if (sum(z$V6=='FDR 1 0p Light' & z$V4=='Entry')) {
# 	# 	rts <-	subset(z, grepl('Port',V6))$time_stamp - subset(z, V6=='FDR 1 2p Light' & V4=='Entry')$time_stamp
# 	# }
# 	# EXTRACT first port entry
# 	rt_position <- min(which(rts>0))
# 	# if length(subset(z, grepl('Port',V6))[1,]==0) {first_port_entry<-0}
# 	first_port_entry <- as.character(subset(z, grepl('Port',V6))$V6[rt_position])
# 	rts <- rts[which(rts>0)] # exclude the early pokes >> animal collecting the previous reward.
# 	if (length(rts)>1) {rts<-rts[1]}
# 	if (length(rts)==0) {rts<-0} # assign 0 to trials with port pickup error.
# 	if (length(first_port_entry)==0) {first_port_entry<-0}
# 	return(c(tp, reward_port, first_port_entry, rts, pellet, pellet_test))
# }, .progress='text')
# colnames(all_data)[c(4:10)] <- c('trial_nr','tp', 'reward_port', 'first_port_entry', 'rt', 'pellet', 'pellet_test')
# all_data$first_port_entry <- ifelse(all_data$first_port_entry=='Port1', 1, ifelse(all_data$first_port_entry=='Port2', 2, 0))
# # combine pellet from 'trianing' and 'test' trials
# #tmp<-all_data
# all_data$pellet <-ifelse(all_data$reward_port==31 & all_data$pellet_test==1, 1,
# 												ifelse(all_data$reward_port==31 & all_data$pellet_test==0, 0,
# 													ifelse(all_data$reward_port==32 & all_data$pellet_test==2, 2,
# 														ifelse(all_data$reward_port==32 & all_data$pellet_test==0, 0, all_data$pellet))))
# write.table(all_data, '/Users/tadeusz/Data/metacog_rat_orsay/data_frame/data_training_step_5.csv')
