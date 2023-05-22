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
  dir <- paste('/Users/PC/Dropbox/eloras_rats/data/step6/S', toString(ss),'/', sep='')
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

write.table(dat, '/Users/PC/Dropbox/eloras_rats/data/deneme/dat_step6expp.csv')


#=======================================================================================
# GET the trial data + UTILS & functions
#=======================================================================================

dat <- read.table('/Users/PC/Dropbox/eloras_rats/data/deneme/dat_step6expp.csv')

all_data <- ddply(dat, .(rat_nr, session, trial_index, group_width, group_dir), function(z){ #int_offset
	#z <- subset(dat, rat_nr==10 & session==3 & trial_index==15)
	tmp1 <- subset(z, V5=='StartPoke' & V3=='Entry')
	tmp2 <- subset(z, V5=='StartPort' & V3=='Input')
	

	start_light = tmp1$time_stamp[1] # nose poke start
	init_entry = tmp2$time_stamp[1] # initial entry in start port

	cond=0
	if (dim(subset(z, V5=='Duration_1s' & V3=='Exit'))[1]>0) {cond=1 # find out group_condition
	} else if (dim(subset(z, V5=='Duration_2s' & V3=='Exit'))[1]>0) {cond=2
	} else if (dim(subset(z, V5=='Duration_4s' & V3=='Exit'))[1]>0) {cond=4}


	int_offset = subset(z, grepl('Duration_', V5) & V3=='Exit')$time_stamp[1]
	first_check = subset(z, grepl('NosePoke', V5) & time_stamp>int_offset)$V4[1] # check first state

	first_check_lat = subset(z, grepl('NosePoke', V5) & time_stamp>int_offset)$time_stamp[1]
	sumNP1 = sum(z$V5 == 'NosePoke1' & z$V3 == 'Input')
	sumNP2 = sum(z$V5 == 'NosePoke2' & z$V3 == 'Input')
	sumNP3 = sum(z$V5 == 'NosePoke3' & z$V3 == 'Input')
	sumNP4 = sum(z$V5 == 'NosePoke4' & z$V3 == 'Input')
	sumNP5 = sum(z$V5 == 'NosePoke5' & z$V3 == 'Input')
	time_NP1 <- z$V2[z$V5 == "NosePoke1" & z$V3 == "Input" & z$V2 > 0]
	time_NP2 <- z$V2[z$V5 == "NosePoke2" & z$V3 == "Input" & z$V2 > 0]
	time_NP3 <- z$V2[z$V5 == "NosePoke3" & z$V3 == "Input" & z$V2 > 0]
	time_NP4 <- z$V2[z$V5 == "NosePoke4" & z$V3 == "Input" & z$V2 > 0]
	time_NP5 <- z$V2[z$V5 == "NosePoke5" & z$V3 == "Input" & z$V2 > 0]
	start_lat <- z$V2[max(which(grepl("Duration_", z$V5)& z$V3 == "Entry"))]
	
	
	lat_NP1 <- min(time_NP1[time_NP1 > start_lat]) - start_lat
	lat_NP2 <- min(time_NP2[time_NP2 > start_lat]) - start_lat
	lat_NP3 <- min(time_NP3[time_NP3 > start_lat]) - start_lat
	lat_NP4 <- min(time_NP4[time_NP4 > start_lat]) - start_lat
	lat_NP5 <- min(time_NP5[time_NP5 > start_lat]) - start_lat
	
	first_poke <- NA
	
	if (!any(is.na(c(lat_NP1, lat_NP2, lat_NP3, lat_NP4, lat_NP5)))) {
	  if (lat_NP1 == min(lat_NP1, lat_NP2, lat_NP3, lat_NP4, lat_NP5)) {
	    first_poke <- 1
	  } else if (lat_NP2 == min(lat_NP1, lat_NP2, lat_NP3, lat_NP4, lat_NP5)) {
	    first_poke <- 2
	  } else if (lat_NP3 == min(lat_NP1, lat_NP2, lat_NP3, lat_NP4, lat_NP5)) {
	    first_poke <- 3
	  } else if (lat_NP4 == min(lat_NP1, lat_NP2, lat_NP3, lat_NP4, lat_NP5)) {
	    first_poke <- 4
	  } else if (lat_NP5 == min(lat_NP1, lat_NP2, lat_NP3, lat_NP4, lat_NP5)) {
	    first_poke <- 5
	  } }
	  
	  
	
	
	return(c(cond, start_light, init_entry, int_offset, first_check, first_check_lat, sumNP1, sumNP2, sumNP3, sumNP4, sumNP5, lat_NP1, lat_NP2, lat_NP3, lat_NP4, lat_NP5, first_poke))
}, .progress='text')

colnames(all_data)[c(6:22)] <- c('cond','start_light', 'init_entry', 'int_offset', 'first_check', 'first_check_lat', 'sumNP1','sumNP2','sumNP3','sumNP4','sumNP5','lat_NP1','lat_NP2','lat_NP3','lat_NP4','lat_NP5', 'first_poke')
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


write.table(all_data, 'C:/Users/PC/Dropbox/eloras_rats/data/deneme/all_data_step6exp.csv')




# Load the required packages
library(tidyr)
library(stringr)

# Read in the CSV file
df <- read.csv("/Users/PC/Dropbox/eloras_rats/data/deneme/all_data_step6exp.csv", header = FALSE, col.names = "data")

# Split the single column into multiple columns
df <- separate(df, data, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16","col17","col18","col19","col20","col21","col22","col23","col24"), sep = " ")



# Convert columns to their appropriate data types
df$col1 <- as.integer(df$col1)
df$col2 <- as.character(df$col2)
df$col3 <- as.integer(df$col3)
df$col4 <- as.integer(df$col4)
df$col5 <- as.character(df$col5)
df$col6 <- as.character(df$col6)
df$col7 <- as.integer(df$col7)
df$col8 <- as.numeric(df$col8)
df$col9 <- as.numeric(df$col9)
df$col10 <- as.numeric(df$col10)
df$col11 <- as.numeric(df$col11)
df$col12 <- as.numeric(df$col12)
df$col13 <- as.numeric(df$col13)
df$col14 <- as.numeric(df$col14)
df$col15 <- as.numeric(df$col15)
df$col16 <- as.numeric(df$col16)
df$col17 <- as.numeric(df$col17)
df$col18 <- as.numeric(df$col18)
df$col19 <- as.numeric(df$col19)
df$col20 <- as.numeric(df$col20)
df$col21 <- as.numeric(df$col21)
df$col22 <- as.numeric(df$col22)
df$col23 <- as.numeric(df$col23)
df$col24 <- as.numeric(df$col24)

# View the resulting data frame
write.csv(df, file = "dnmall_data_step6expcolumned.csv", row.names = FALSE)






# read in the data from the CSV file and exclude lat_NP values greater than 20
all_data <- subset(read.csv('C:/Users/PC/Dropbox/eloras_rats/data/deneme/dnmall_data_step6expcolumned.csv'))

# create four subsets based on the given conditions
subset1 <- all_data[all_data$cond == 4 & all_data$group_width == "NARO" & all_data$group_dir == "LM", ]
subset2 <- all_data[all_data$cond == 4 & all_data$group_width == "NARO" & all_data$group_dir == "ML", ]
subset3 <- all_data[all_data$cond == 4 & all_data$group_width == "WIDE" & all_data$group_dir == "LM", ]
subset4 <- all_data[all_data$cond == 4 & all_data$group_width == "WIDE" & all_data$group_dir == "ML", ]

# create SE and ME columns for each subset
subset1$SE <- ifelse(is.infinite(subset1$lat_NP3), NA, subset1$lat_NP3)
subset1$ME <- ifelse(is.infinite(subset1$lat_NP4), NA, subset1$lat_NP4)
subset2$SE <- ifelse(is.infinite(subset2$lat_NP3), NA, subset2$lat_NP3)
subset2$ME <- ifelse(is.infinite(subset2$lat_NP2), NA, subset2$lat_NP2)
subset3$SE <- ifelse(is.infinite(subset3$lat_NP3), NA, subset3$lat_NP3)
subset3$ME <- ifelse(is.infinite(subset3$lat_NP5), NA, subset3$lat_NP5)
subset4$SE <- ifelse(is.infinite(subset4$lat_NP3), NA, subset4$lat_NP3)
subset4$ME <- ifelse(is.infinite(subset4$lat_NP1), NA, subset4$lat_NP1)

# For subset 1
subset1$SE <- ifelse(subset1$first_poke == 3, 1, 0)
subset1$ME <- ifelse(subset1$first_poke == 4, 1, 0)

# For subset 2
subset2$SE <- ifelse(subset2$first_poke == 3, 1, ifelse(subset2$first_poke == 2, 0, NA))
subset2$ME <- ifelse(subset2$first_poke == 2, 1, 0)

# For subset 3
subset3$SE <- ifelse(subset3$first_poke == 3, 1, 0)
subset3$ME <- ifelse(subset3$first_poke == 5, 1, 0)

# For subset 4
subset4$SE <- ifelse(subset4$first_poke == 3, 1, 0)
subset4$ME <- ifelse(subset4$first_poke == 1, 1, 0)
# combine the subsets back into one data frame
all_data <- rbind(subset1, subset2, subset3, subset4)

# write the modified data frame to a new CSV file with lat_NP columns
write.csv(all_data, 'C:/Users/PC/Dropbox/eloras_rats/data/deneme/all_data_step6exp_filtered.csv', row.names = FALSE)




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
