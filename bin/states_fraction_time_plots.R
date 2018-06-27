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
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. June 2018                  ###
#############################################################################
### Plotting fraction of time of each state from data segmented using     ###
### chromHMM.  
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

### Execution example
## Rscript plot_HMM_segmentation.R --path_bed_files="path_to_bed_files"
library("ggplot2")
library ("plyr") # ddply
library ("dplyr") # ddply
# library("reshape2")
# library ("gtools")
# library("tidyr")

# Loading params plot:
source("https://raw.githubusercontent.com/cbcrg/mwm/master/lib/R/plot_param_public.R")

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
      states_fraction_time_plots.R
      
      Arguments:
      --path_to_tbl=path_to_tbl         - character
      --path_to_tbl_col=path_to_tbl_col - character
      --image_format=image_format       - character
      --help                            - print this text
      
      Example:      
      ./states_fraction_time_plots.R --path_to_tbl=\"path_to_tbl\" --path_to_tbl_col=\"path_to_tbl_col\" --image_format=\"image_format\" \n")
  
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

# All arguments are mandatory
{
  if (is.null (argsL$path_to_tbl))
  {
    stop ("[FATAL]: path to tbl containing chromHMM segmentation not provided", stderr())
  }
  else
  {
    path_to_tbl <- argsL$path_to_tbl
  }
}

# tbl of colors
{
  if (is.null (argsL$path_to_tbl_col))
  {
    stop ("[FATAL]: path to tbl containing color of HMM states not provided", stderr())
  }
  else
  {
    path_to_tbl_col <- argsL$path_to_tbl_col
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

df_fraction_t_by_phase <- read.csv (path_to_tbl, header=F, sep="\t")
colnames(df_fraction_t_by_phase) <- c("state_n", "start", "end", "phase", "group", "mouse_id", 
                                      "start2", "end2", "color","counts", "fraction", "total", 
                                      "percentage", "day_phase")

color_state <- read.csv(path_to_tbl_col, header=F, sep="\t")
colnames(color_state) <- c("state_n", "color")

color_state$color_code <- unlist(lapply(as.vector(color_state$color), function (x) sapply(strsplit(x, ","), function(y)
  rgb(y[1], y[2], y[3], maxColorValue=255))))

summary_fract <- ddply(df_fraction_t_by_phase, .(state_n, group, phase, day_phase), summarise, mean=mean(percentage))
summary_fract <- reshape(summary_fract, idvar=c("state_n", "group",  "day_phase"), timevar="phase", direction="wide")
colnames(summary_fract)[4:5] <- c("Development", "Habituation")

ctrl_df_light <- subset (summary_fract, group=="Ctrl" & day_phase=="light")
ctrl_df_light$cum_sum <- cumsum(ctrl_df_light$Habituation)
ctrl_df_dark <- subset (summary_fract, group=="Ctrl" & day_phase=="dark")
ctrl_df_dark$cum_sum <- cumsum(ctrl_df_dark$Habituation)
hf_df_light <- subset (summary_fract, group=="HF" & day_phase=="light")
hf_df_light$cum_sum <- cumsum(hf_df_light$Habituation)
hf_df_dark <- subset (summary_fract, group=="HF" & day_phase=="dark")
hf_df_dark$cum_sum <- cumsum(hf_df_dark$Habituation)

spie.data <- rbind (ctrl_df_light, ctrl_df_dark, hf_df_light, hf_df_dark)
spie.data$state_n <- as.character(spie.data$state_n)
levels(spie.data$day_phase) <- c("Dark phase", "Light phase")
spie.data$day_phase <- factor(spie.data$day_phase, levels = c("Light phase", "Dark phase"))

bars <-ggplot(spie.data, aes(x = cum_sum - 0.5 * Habituation, fill = state_n), labels=state_n) +
  geom_bar(aes(width = Habituation, y = 1),
           color = NA, stat = "identity", alpha = 0.2) +
  geom_bar(aes(width = Habituation, y = Development / Habituation),
           color = "grey10", size = 0.1, stat = "identity") +
  annotate("segment", linetype = 3, size = 0.5, lineend = "round",
           x = -Inf, xend = Inf, y = 1.0, yend = 1.0) +
  scale_fill_manual(values = color_state$color_code,
                    labels=c("Long meals", "Regular meals (2/3 min)", "Inactive", "Short meals")) +
  scale_y_sqrt(limits = c(0, 8)) +
  # scale_x_discrete (labels = state_n,
                      # scale_x_continuous (labels = "",
                      # breaks = cum_sum) +
  # geom_text(data=spie.data, aes(label=p_value, y = Development / Habituation), vjust=-0.25) +
  facet_grid(group~day_phase)
# print(bars)

size_strips <- 12

spie <- bars + coord_polar (theta = "x") + #scale_y_log10()+
  labs (x = "State", 
        y = "",
        title = "Distribution of HMM states during development\nphase relative to habituation phase\n", 
        fill = "State") + 
  theme_bw() +
  # theme (axis.text.x = element_text(angle = 0, size = 8)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_text(size=size_strips, face="bold"), 
        strip.text.y = element_text(size=size_strips, face="bold"),
        panel.border = element_rect(colour = "black"))
  
print (spie)

plot_width <- 14
plot_height <- 8

name_out <- paste ("states_fraction", ".", image_format, sep="")
ggsave (file=name_out, width = plot_width, height=plot_height)

## Significance
# t_test_result <- df_fraction_t_by_phase %>% 
#   group_by(group, day_phase, state_n) %>%
#   summarise(pval = t.test(percentage ~ phase, var.equal = TRUE)$p.value)
# 
# spie.data$p_value <- t_test_result$pval
# signif <- symnum(t_test_result$pval, corr = FALSE, na = FALSE, 
#                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
# spie.data$p_value <- signif
