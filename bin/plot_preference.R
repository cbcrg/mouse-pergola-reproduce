#!/usr/bin/env Rscript

#  Copyright (c) 2014-2017, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2017, Jose Espinosa-Carrasco and the respective authors.
#
#  This file is part of Pergola.
#
#  Pergola is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Pergola is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. June 2016                  ###
#############################################################################
### Mean values for each group of mice                                    ###
### Using bed files with bouts intersected with experimental              ###
### phases or other phases                                                ###
#############################################################################

##Getting HOME directory
# home <- Sys.getenv("HOME") #del

##Loading libraries
library ("ggplot2")
# library ("plotrix") #std.error
# library('extrafont')
library ('gtools') 
library(dplyr)

#####################
### VARIABLES
#Reading arguments
args <- commandArgs (TRUE) #if not it doesn't start to count correctly

## Default setting when no arguments passed
if ( length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      starting_regions_file_vs_24h
      
      Arguments:
      --stat=someValue        - character, stat to analyze (sum, mean, ...)
      --path2files=someValue - character, path to read files
      --help                 - print this text
      
      Example:
      ./starting_regions_file_vs_24h.R --stat=\"sum\" \n")
  
  q (save="no")
}

# Use to parse arguments beginning by --
parseArgs <- function(x) 
{
  strsplit (sub ("^--", "", x), "=")
}

#Parsing arguments
argsDF <- as.data.frame (do.call("rbind", parseArgs(args)))
argsL <- as.list (as.character(argsDF$V2))
names (argsL) <- argsDF$V1

# stat is mandatory
{
  if (is.null (argsL$stat)) 
  {
    stop ("[FATAL]: Stat parameter is mandatory")
  }
  else
  {
    stat <- argsL$stat
  }
}

# path to files
{
  if (is.null (argsL$path2files)) 
  {
    stop ("[FATAL]: Path to files is mandatory")
  }
  else
  {
    path2files <- argsL$path2files
  }
}

## Loading parameters for plots
source("https://gist.githubusercontent.com/JoseEspinosa/307430167b05e408ac071b8724701abf/raw/06b26f764953ceea53d334190d7b736308e5141d/param_plot_publication.R")

# path2files <- "/Users/jespinosa/git/pergola/examples/CB1_mice_tt/results_by_day/feeding_by_phases/sum/" #del
pwd <- getwd()
setwd(path2files)

files <- list.files(pattern=paste("tr_.*.bed$", sep=""))
# stat<-"sum"
data.frame_bed <- NULL

for (bed_file in files) {

  info = file.info (bed_file)
  if (info$size == 0) { next }
  df <- read.csv(bed_file, header=F, sep="\t")
  phenotype <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[1])
  genotype <- unlist(strsplit(phenotype, split="_",fixed=T))[1]
  mouse <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[2])
  data_type <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[3])
  data_type <- gsub ("food_fat", "HF", data_type)
  data_type <- gsub ("food_sc", "SC", data_type)
  phase <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[4])
  exp_phase <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[5])
  df$phenotype <- phenotype
  df$phenotype <- factor(df$phenotype, levels=c("wt_saline", "wt_nicotine", "KO_cb1_saline", "KO_cb1_nicotine"),
                         labels=c("wt_saline", "wt_nicotine", "KO_cb1_saline", "KO_cb1_nicotine"))
#   df$phenotype <- factor(df$phenotype, levels=c("WT Saline", "WT Nicotine", "KO Saline", "KO Nicotine"),
#                                                 labels=c("WT Saline", "WT Nicotine", "KO Saline", "KO Nicotine"))
                                                  
  df$mouse <- mouse
  df$genotype <- genotype
  df$genotype <- factor(df$genotype, levels=c("wt", "KO"),
                       labels=c("wt", "KO"))
#   df$genotype <- factor(df$genotype, levels=c("wt", "KO"),
#                       labels=c("wt", "KO"))
  df$mouse <- mouse
  df$data_type <- data_type
  df$phase <- phase
  df$exp_phase <- gsub("_", " ", exp_phase)
  df$group2plot <- paste (phase, data_type)
  data.frame_bed <- rbind(data.frame_bed, df)
}
setwd(pwd)

## animal 29 has not data for nicotine treatment
## 3799 chr1  0  1 no_hits  0  0     wt_saline    29        SC   day  Nicotine treatment     day SC
to_replace <- subset(data.frame_bed, mouse=='1' & exp_phase=="Nicotine treatment" & data_type=="SC")
to_replace$V6 <- 0
data.frame_bed <- rbind(to_replace, data.frame_bed)

#######
# In the analysis of the original paper they substitute empty values for an animal
# for the medians, but not only for these days but for the whole experimental
# period during which the empty values were occurring
# This part of the code will be use to preprocess the data and use the same approach
# in order to reproduce the plots they published

nzmean <- function(x) {
  zvals <- x==0
  if (all(zvals)) 0 else mean(x[!zvals])
  # if (all(zvals)) 0 else median(x[!zvals])
}

means_no_zeros <- with (data.frame_bed, aggregate (cbind (V6), list (genotype=genotype, data_type=data_type,
                        exp_phase=exp_phase, day=V4 ),
                        FUN=function (x) c (mean=nzmean(x))))

for (i in c(1:length (data.frame_bed$exp_phase))) {
  if (data.frame_bed [i,"V6"] == 0) {
    # If I have modified one day the rest are not anymore 0
    # I have to change by the median value in the rest of this day
    data.frame_bed [i,"V6"] <- subset(means_no_zeros, genotype==data.frame_bed [i,"genotype"] &
                                      exp_phase==data.frame_bed [i,"exp_phase"] &
                                      data_type==data.frame_bed [i,"data_type"] &
                                      day==data.frame_bed [i,"V4"], select=c("V6"))
  }
}

## delete the first 4 days of the nicotine tt -> injection
data.frame_bed <- subset (data.frame_bed, V4!="day_15" & V4!="day_16" & V4!="day_17" & V4!="day_18")

## color blind friendly palette
{
  if ('HF' %in% data.frame_bed$data_type & 'SC' %in% data.frame_bed$data_type) {
    cbb_palette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    bouts_title <- "Feeding bouts"
    units <- 'g'
    filter_pref <- 'HF'
  }
  else if ('water' %in% data.frame_bed$data_type & 'saccharin' %in% data.frame_bed$data_type) {
    cbb_palette <- c("#0072B2", "#D55E00", "#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
    bouts_title <- "Drinking bouts"
    units <- 'mL'
    filter_pref <- 'saccharin'
  }
  else {
    cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    bouts_title <- "Number of bouts"
    units <- 'no units'
  }
}
size_titles <- 20
size_axis <- 18
size_axis_ticks_x <- 14
size_axis_ticks_y <- 14
size_strip_txt <- 14
plot_width <- 12
plot_height <- 10
font <- "Arial"

name_file <- "plot"

# plot_title <- switch (stat, count=bouts_title, mean="Mean intake per feeding bout", median='Median intake per feeding bout', sum='Accumulated intake',
#                       max='Maximun intake', '')
plot_title <- "Preference"
## Combined phases as in the original paper
tbl_stat_mean_comb_ph <- with (data.frame_bed, aggregate (cbind (V6), list (phenotype=phenotype, data_type=data_type,
                                                                    exp_phase=exp_phase),
                                                  FUN=function (x) c (mean=mean(x))))

tbl_stat_mean_comb_ph$mean <- tbl_stat_mean_comb_ph$V6

######################
# preference by animal
#preference <- mean_intake HF  / (mean_intake SC + mean_intake HF) They do animal by animal and file by file

# same thing but only with basal
data.frame_bed_basal <- subset (data.frame_bed, exp_phase=="Basal")

{
  if (stat == 'count' | stat == 'sum'){

    data.frame_bed$value <- data.frame_bed$V6
    preference <- data.frame_bed %>%
      # they do in functin of the last one (data types)
      group_by(mouse, phenotype, exp_phase, phase,  V2, V3, data_type) %>%
      summarise(mean=mean(value)) %>%
      mutate(pref=mean/sum(mean)*100)

    preference <- as.data.frame(preference)
    preference <- preference[!is.na(preference$pref),]

    preference_mean_comb_ph <- with (preference, aggregate (cbind (pref), list (phenotype=phenotype, data_type=data_type,
                                                                                exp_phase=exp_phase),
                                                              FUN=function (x) c (mean=mean(x))))

    preference_mean_comb_ph$mean <- preference_mean_comb_ph$pref

    ## basal period
    data.frame_bed_basal$value <- data.frame_bed_basal$V6
    preference_basal <- data.frame_bed_basal %>%
      # they do in function of the last one (data types)
      group_by(mouse, genotype, exp_phase, phase,  V2, V3, data_type) %>%
      summarise(mean=mean(value)) %>%
      mutate(pref=mean/sum(mean)*100)
    
    preference_basal <- as.data.frame(preference_basal)
    preference_basal <- preference_basal[!is.na(preference_basal$pref),]
    preference_mean_comb_ph_basal <- with (preference_basal, aggregate (cbind (pref), list (genotype=genotype, data_type=data_type,
                                                                                exp_phase=exp_phase),
                                                            FUN=function (x) c (mean=mean(x))))

    preference_mean_comb_ph_basal$mean <- preference_mean_comb_ph_basal$pref

    name_file <- paste ("preference", ".", "png", sep="")

    ## Join the two tables
    preference_mean_comb_ph <- subset (preference_mean_comb_ph, exp_phase!="Basal")
    colnames(preference_mean_comb_ph_basal)[1] <- "phenotype"
    preference_mean_comb_ph_plot <- rbind(preference_mean_comb_ph_basal, preference_mean_comb_ph)

    ggplot(data=preference_mean_comb_ph_plot, aes(x=phenotype, y=pref, fill=data_type)) +
      geom_bar(stat="identity", position="fill") +
      scale_fill_manual(values = c(cbb_palette[2], cbb_palette[1])) +
      scale_y_continuous (breaks=seq(0, 1,0.2)) +
      labs (title = paste(plot_title, "\n", sep="")) +
      labs (y = paste(paste ("Percentage", "\n", sep="")), x="\n") +
#       theme (plot.title = element_text(family=font, size=size_titles)) +
#       theme (axis.title.x = element_text(family=font, size=size_axis)) +
#       theme (axis.title.y = element_text(family=font, size=size_axis)) +
#       theme (axis.text.x = element_text(family=font, size=size_axis_ticks_x, angle=90)) +
#       theme (axis.text.y = element_text(family=font, size=size_axis_ticks_y)) +
#       theme (axis.text.x = element_text(family=font, angle=90, vjust=0.4,hjust=1)) +
      theme (plot.title = element_text(size=size_titles)) +
      theme (axis.title.x = element_text(size=size_axis)) +
      theme (axis.title.y = element_text(size=size_axis)) +
      theme (axis.text.x = element_text(size=size_axis_ticks_x, angle=90)) +
      theme (axis.text.y = element_text(size=size_axis_ticks_y)) +
      theme (axis.text.x = element_text(angle=90, vjust=0.4,hjust=1)) +  
      facet_wrap(~exp_phase, ncol=3, scale="free") +
      scale_x_discrete(labels=c("WT", "KO", "WT Saline", "WT Nicotine", "KO Saline", "KO Nicotine", 
                                "WT Saline", "WT Nicotine", "KO Saline", "KO Nicotine")) +
      theme(strip.background = element_rect(fill="white")) +
      theme(strip.text.x = element_text(size = size_axis_ticks_x), legend.title=element_blank())
    
    ggsave (file=name_file, width=plot_width, height=plot_height, dpi=300)

  }
}

