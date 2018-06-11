#!/usr/bin/env Rscript

#  Copyright (c) 2014-2018, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2018, Jose Espinosa-Carrasco and the respective authors.
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
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. Jan 2018                   ###
#############################################################################
### Heatmap comparing fedding behavior                                    ###
#############################################################################

# library (gtools) #foldchange
# library (plyr) 
library (ggplot2)

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
      plot_em_trans_by_win.R
      
      Arguments:
      --path2ems=someValue        - character, path to read emissions table
      --path2trans=someValue      - character, path to read transitions table
      --image_format=image_format - character
      --help                      - print this text
      
      Example:
      ./.R --path2ems=\"path_to_emission_tbl\" --path2trans=\"path_to_transitions_tbl\" --image_format=\"image_format\" \n")
  
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

# path to emissions
{
  if (is.null (argsL$path2ems))
  {
    stop ("[FATAL]: Path to emissions table is mandatory")
  }
  else
  {
    path2ems <- argsL$path2ems
  }
}

# path to transitions
{
  if (is.null (argsL$path2trans))
  {
    stop ("[FATAL]: Path to transition table is mandatory")
  }
  else
  {
    path2trans <- argsL$path2trans
  }
}

# plot image format
{
  if (is.null (argsL$image_format))
  {
    image_format <- "tiff"
    warning ("[Warning]: format for plots not provided, default tiff")
  }
  else
  {
    image_format <- argsL$image_format
  }
}

source("https://gist.githubusercontent.com/JoseEspinosa/307430167b05e408ac071b8724701abf/raw/06b26f764953ceea53d334190d7b736308e5141d/param_plot_publication.R")

pwd <- getwd()

## color blind friendly palette
cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# path2ems <- "/Users/jespinosa/Desktop/tbls_HMM/emission_prob.tbl"
# path2ems <- "/Users/jespinosa/git/cbcrg/mouse-pergola-reproduce/results/chromHMM/emission_prob_win.tbl"
# path2ems <- "/Users/jespinosa/git/cbcrg/mouse-pergola-reproduce/results/chromHMM/emission_prob.tbl"
data.frame_prob_emission <- read.csv(path2ems, header=F, sep="\t")

colnames(data.frame_prob_emission) <- c("day", "group", "prob_mode", "State", "Emission_n", "Emission", "binary", "prob")
data.frame_prob_emission_filt <- data.frame_prob_emission [which(data.frame_prob_emission$binary==1), ]

### CTRL
## 2 to 1 and 1 to 2
days_crossed_ctrl <- data.frame_prob_emission_filt$day [ which (data.frame_prob_emission_filt$State == 2 & 
                                                                  data.frame_prob_emission_filt$group == "ctrl" &
                                                                  data.frame_prob_emission_filt$Emission == "120_food" &
                                                                  data.frame_prob_emission_filt$prob > 0.2)] 

data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_ctrl &
                                               data.frame_prob_emission_filt$State == 2 & 
                                               data.frame_prob_emission_filt$group == "ctrl")] <- "s1"
data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_ctrl &
                                               data.frame_prob_emission_filt$State == 1 & 
                                               data.frame_prob_emission_filt$group == "ctrl")] <- "s2"
data.frame_prob_emission_filt$State <- gsub("s", "", data.frame_prob_emission_filt$State)

## 3 to 1 and 1 to 3
days_crossed_ctrl_1 <- data.frame_prob_emission_filt$day [ which (data.frame_prob_emission_filt$State == 3 & 
                                                                    data.frame_prob_emission_filt$group == "ctrl" &
                                                                    data.frame_prob_emission_filt$Emission == "120_food" &
                                                                    data.frame_prob_emission_filt$prob > 0.2)] 

data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_ctrl_1 &
                                               data.frame_prob_emission_filt$State == 3 & 
                                               data.frame_prob_emission_filt$group == "ctrl")] <- "s1"
data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_ctrl_1 &
                                               data.frame_prob_emission_filt$State == 1 & 
                                               data.frame_prob_emission_filt$group == "ctrl")] <- "s3"
data.frame_prob_emission_filt$State <- gsub("s", "", data.frame_prob_emission_filt$State)

## 3 to 2 and 2 to 3
days_crossed_ctrl_2 <- data.frame_prob_emission_filt$day [ which (data.frame_prob_emission_filt$State == 3 &
                                                                    data.frame_prob_emission_filt$group == "ctrl" &
                                                                    data.frame_prob_emission_filt$Emission == "30_120_food" &
                                                                    data.frame_prob_emission_filt$prob < 0.001) ]

data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_ctrl_2 &
                                               data.frame_prob_emission_filt$State == 3 &
                                               data.frame_prob_emission_filt$group == "ctrl")] <- "s2"
data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_ctrl_2 &
                                               data.frame_prob_emission_filt$State == 2 &
                                               data.frame_prob_emission_filt$group == "ctrl")] <- "s3"
data.frame_prob_emission_filt$State <- gsub("s", "", data.frame_prob_emission_filt$State)

### HF
days_crossed_hf <- data.frame_prob_emission_filt$day [ which (data.frame_prob_emission_filt$State == 1 &
                                                                data.frame_prob_emission_filt$group == "hf" &
                                                                data.frame_prob_emission_filt$Emission == "120_food" &
                                                                data.frame_prob_emission_filt$prob < 0.5) ]
## 3 to 1 and 2 to 3
data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_hf &
                                               data.frame_prob_emission_filt$State == 3 &
                                               data.frame_prob_emission_filt$group == "hf")] <- "s1"
data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_hf &
                                               data.frame_prob_emission_filt$State == 1 &
                                               data.frame_prob_emission_filt$group == "hf")] <- "s3"
data.frame_prob_emission_filt$State <- gsub("s", "", data.frame_prob_emission_filt$State)

data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_hf &
                                               data.frame_prob_emission_filt$State == 2 &
                                               data.frame_prob_emission_filt$group == "hf")] <- "s3"
data.frame_prob_emission_filt$State [ which (data.frame_prob_emission_filt$day %in%  days_crossed_hf &
                                               data.frame_prob_emission_filt$State == 3 &
                                               data.frame_prob_emission_filt$group == "hf")] <- "s2"
data.frame_prob_emission_filt$State <- gsub("s", "", data.frame_prob_emission_filt$State)

data.frame_prob_emission_filt$emission_translated <- gsub ("0_30_food", "Short meals", data.frame_prob_emission_filt$Emission)
data.frame_prob_emission_filt$emission_translated <- gsub ("30_120_food", "Mid meals", data.frame_prob_emission_filt$emission_translated)
data.frame_prob_emission_filt$emission_translated <- gsub ("120_food", "Long meals", data.frame_prob_emission_filt$emission_translated)
data.frame_prob_emission_filt$emission_translated <- gsub ("30_food", "Long meals", data.frame_prob_emission_filt$emission_translated)

data.frame_prob_emission_filt$state_emission <- paste(data.frame_prob_emission_filt$State, data.frame_prob_emission_filt$emission_translated,  sep="_")

# ggplot(data.frame_prob_emission_filt, aes(x=day, y=prob, group=state_emission,  colour=state_emission)) +
#   geom_line() +
#   facet_grid(group ~ state_emission) +
#   labs(title = "States emissions \n(1 week window with 1 day step)", x = "\nDays", y="Probability\n") +
#   theme(plot.title = element_text(hjust = 0.5)) 
xmin <- 0
xmax <- 62
xmin <- 0
xmax <- 8

ymin <- 0
ymax <- 1

ggplot(data.frame_prob_emission_filt, aes(x=day, y=prob, group=group,  colour=group)) + 
  geom_line() +
  # facet_grid(. ~ state_emission) +
  facet_grid(State ~ emission_translated) +
  labs(title = "States emissions (1 week window with 1 day step)\n", x = "\nDays", y="Probability\n") +
  scale_x_continuous (limits=c(xmin, xmax)) +
  scale_y_continuous (limits=c(ymin, ymax)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = cbb_palette, name='', labels = c("Ctrl", "HF")) +
  theme (legend.key.size = unit(0.8, "cm"), legend.title = element_blank()) # +
  # geom_vline(xintercept = 6, linetype="dashed") +
  # geom_vline(xintercept = 12, linetype="dashed")
# 4*7 --> 28
# 6*7 --> 42
plot_width <- 13
plot_height <- 8

# name_out <- "/Users/jespinosa/Desktop/emission.png"
name_out <- paste ("emission", ".", image_format, sep="")
ggsave (file=name_out, width = plot_width, height=plot_height)

## TRANSITIONS
# path2trans <- "/Users/jespinosa/Desktop/tbls_HMM/transition_prob.tbl"

data.frame_prob_trans <- read.csv(path2trans, header=F, sep="\t")
head (data.frame_prob_trans)

# path2trans <- "/Users/jespinosa/git/cbcrg/mouse-pergola-reproduce/results/chromHMM/transition_prob.tbl"
data.frame_prob_trans_2 <- read.csv(path2trans, header=F, sep="\t")
summary(data.frame_prob_trans_2$V6)
summary(data.frame_prob_trans$V6)

colnames(data.frame_prob_trans) <- c("day", "group", "prob_mode", "State_ini", "State_end", "prob")

### CTRL
## 2 to 1 and 1 to 2
## days_crossed_ctrl
data.frame_prob_trans$State_ini [ which (data.frame_prob_trans$day %in% days_crossed_ctrl &
                                           data.frame_prob_trans$group == "ctrl" &
                                           data.frame_prob_trans$State_ini == 1 ) ] <- "s2"  
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_ctrl &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_end == 1 )] <- "s2"  

data.frame_prob_trans$State_ini [which (data.frame_prob_trans$day %in% days_crossed_ctrl &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_ini == 2 )] <- "s1"  
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_ctrl &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_end == 2 )] <- "s1"  

data.frame_prob_trans$State_ini <- gsub("s", "", data.frame_prob_trans$State_ini)
data.frame_prob_trans$State_end <- gsub("s", "", data.frame_prob_trans$State_end)

## 3 to 1 and 1 to 3
# days_crossed_ctrl_1 
data.frame_prob_trans$State_ini [ which (data.frame_prob_trans$day %in% days_crossed_ctrl_1 &
                                           data.frame_prob_trans$group == "ctrl" &
                                           data.frame_prob_trans$State_ini == 3 ) ] <- "s1"
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_ctrl_1 &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_end == 3 )] <- "s1"

data.frame_prob_trans$State_ini [which (data.frame_prob_trans$day %in% days_crossed_ctrl_1 &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_ini == 1 )] <- "s3"
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_ctrl_1 &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_end == 1 )] <- "s3"

data.frame_prob_trans$State_ini <- gsub("s", "", data.frame_prob_trans$State_ini)
data.frame_prob_trans$State_end <- gsub("s", "", data.frame_prob_trans$State_end)

## 3 to 2 and 2 to 3
# days_crossed_ctrl_2
data.frame_prob_trans$State_ini [ which (data.frame_prob_trans$day %in% days_crossed_ctrl_2 &
                                           data.frame_prob_trans$group == "ctrl" &
                                           data.frame_prob_trans$State_ini == 3 ) ] <- "s2"
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_ctrl_2 &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_end == 3 )] <- "s2"

data.frame_prob_trans$State_ini [which (data.frame_prob_trans$day %in% days_crossed_ctrl_2 &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_ini == 2 )] <- "s3"
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_ctrl_2 &
                                          data.frame_prob_trans$group == "ctrl" &
                                          data.frame_prob_trans$State_end == 2 )] <- "s3"

data.frame_prob_trans$State_ini <- gsub("s", "", data.frame_prob_trans$State_ini)
data.frame_prob_trans$State_end <- gsub("s", "", data.frame_prob_trans$State_end)

data.frame_prob_trans [ which (data.frame_prob_trans$day %in%  days_crossed_hf &
                                 data.frame_prob_trans$group == "hf" &
                                 data.frame_prob_trans$State_ini == 2 ), ]

data.frame_prob_trans$State_ini [ which (data.frame_prob_trans$day %in% days_crossed_hf &
                                           data.frame_prob_trans$group == "hf" &
                                           data.frame_prob_trans$State_ini == 1 ) ] <- "s3"  

data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_hf &
                                          data.frame_prob_trans$group == "hf" &
                                          data.frame_prob_trans$State_end == 1 )] <- "s3"  

data.frame_prob_trans$State_ini [which (data.frame_prob_trans$day %in% days_crossed_hf &
                                          data.frame_prob_trans$group == "hf" &
                                          data.frame_prob_trans$State_ini == 3 )] <- "s1"  
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_hf &
                                          data.frame_prob_trans$group == "hf" &
                                          data.frame_prob_trans$State_end == 3 )] <- "s1"  

data.frame_prob_trans$State_ini <- gsub("s", "", data.frame_prob_trans$State_ini)
data.frame_prob_trans$State_end <- gsub("s", "", data.frame_prob_trans$State_end)

data.frame_prob_trans$State_ini [ which (data.frame_prob_trans$day %in% days_crossed_hf &
                                           data.frame_prob_trans$group == "hf" &
                                           data.frame_prob_trans$State_ini == 2 ) ] <- "s3"  

data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_hf &
                                          data.frame_prob_trans$group == "hf" &
                                          data.frame_prob_trans$State_end == 2 )] <- "s3"  

data.frame_prob_trans$State_ini [which (data.frame_prob_trans$day %in% days_crossed_hf &
                                          data.frame_prob_trans$group == "hf" &
                                          data.frame_prob_trans$State_ini == 3 )] <- "s2"  
data.frame_prob_trans$State_end [which (data.frame_prob_trans$day %in% days_crossed_hf &
                                          data.frame_prob_trans$group == "hf" &
                                          data.frame_prob_trans$State_end == 3 )] <- "s2"  

data.frame_prob_trans$State_ini <- gsub("s", "", data.frame_prob_trans$State_ini)
data.frame_prob_trans$State_end <- gsub("s", "", data.frame_prob_trans$State_end)

data.frame_prob_trans$transition <- paste(data.frame_prob_trans$State_ini, data.frame_prob_trans$State_end,  sep="_")

# ggplot(data.frame_prob_trans, aes(x=day, y=prob, group=transition,  colour=transition)) + 
#   geom_line() +
#   facet_grid(group ~ transition)

xmin <- 0
xmax <- 62
xmin <- 0
xmax <- 8

ymin <- 0
ymax <- 1
  
ggplot(data.frame_prob_trans, aes(x=day, y=prob, group=group,  colour=group)) + 
  geom_line() +
  facet_grid(State_ini ~ State_end) +
  labs(title = "State transitions (1 week window with 1 day step)\n", x = "\nDays", y="Probability\n") +
  # scale_x_continuous (breaks=breaks_v, limits=c(xmin, xmax)) +
  scale_x_continuous (limits=c(xmin, xmax)) +
  scale_y_continuous (limits=c(ymin, ymax)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = cbb_palette, name='', labels = c("Ctrl", "HF")) +
  theme (legend.key.size = unit(0.8, "cm"), legend.title = element_blank()) #+
  # geom_vline(xintercept = 6, linetype="dashed") +
  # geom_vline(xintercept = 42, linetype="dashed")

plot_width <- 13
plot_height <- 8

name_out <- paste ("transition", ".", image_format, sep="")
ggsave (file=name_out, width = plot_width, height=plot_height)

