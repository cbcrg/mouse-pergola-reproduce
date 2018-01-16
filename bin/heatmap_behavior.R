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

library (gtools) #foldchange
library (plyr) 
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
      heatmap_behavior.R

      Arguments:
      --path2files=someValue      - character, path to read files
      --image_format=image_format - character
      --help                      - print this text

      Example:
      ./heatmap_behavior.R --path2files=\"path_to_files\" --image_format=\"image_format\" \n")

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
setwd(path2files)
files <- list.files(pattern=paste("tr_.*.tbl$", sep=""))
data.frame_bed <- NULL

for (bed_file in files) {
  
  info = file.info (bed_file)
  if (info$size == 0) { next }
  df <- read.csv(bed_file, header=F, sep="\t")
  
  df$id <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[1])
  df$group <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[2])
  df$data_type <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[3])
  df$variable <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[4])
  
  data.frame_bed <- rbind(data.frame_bed, df)
}

### Calculate mean by group, variable and week
data.frame_bed.ctrl <- ddply(data.frame_bed [ which (data.frame_bed$group == "ctrl" ), ], .(group, data_type,variable, V4), summarise, mean=mean(V6))
data.frame_bed.hf <- ddply(data.frame_bed [ which (data.frame_bed$group == "hf" ), ], .(group, data_type, variable, V4), summarise, mean=mean(V6))
data.frame_bed.hf$foldChange <- foldchange (data.frame_bed.hf$mean,data.frame_bed.ctrl$mean)

data.frame_bed.hf$period <- as.numeric(gsub ("week_", "", as.character(data.frame_bed.hf$V4)))
colnames(data.frame_bed.hf)<- c("group", "channel", "variable", "week", "value", "foldChange","period")

data.frame_bed.hf$varOrder<- "dummy"
data.frame_bed.hf$varOrder [which (data.frame_bed.hf$variable == "intermeal_time")] <-  "a"
data.frame_bed.hf$varOrder [which (data.frame_bed.hf$variable == "rate")] <-  "b"
data.frame_bed.hf$varOrder [which (data.frame_bed.hf$variable == "duration")] <-  "c"
data.frame_bed.hf$varOrder [which (data.frame_bed.hf$variable == "n_bouts")] <-  "d"
data.frame_bed.hf$varOrder [which (data.frame_bed.hf$variable == "mean")] <-  "e"

data.frame_bed.hf$orderOut<- "dummy"
data.frame_bed.hf$orderOut [which (data.frame_bed.hf$variable == "intermeal_time")] <-  "1"
data.frame_bed.hf$orderOut [which (data.frame_bed.hf$variable == "rate")] <-  "2"
data.frame_bed.hf$orderOut [which (data.frame_bed.hf$variable == "duration")] <-  "3"
data.frame_bed.hf$orderOut [which (data.frame_bed.hf$variable == "n_bouts")] <-  "4"
data.frame_bed.hf$orderOut [which (data.frame_bed.hf$variable == "mean")] <-  "5"

data.frame_bed.hf$variable <- gsub("rate", "Eating rate", data.frame_bed.hf$variable)
data.frame_bed.hf$variable <- gsub("duration", "Average duration", data.frame_bed.hf$variable)
data.frame_bed.hf$variable <- gsub("n_bouts", "Number of meals", data.frame_bed.hf$variable)
data.frame_bed.hf$variable <- gsub("mean", "Average intake", data.frame_bed.hf$variable)
data.frame_bed.hf$variable <- gsub("intermeal_time", "Average intermeal duration", data.frame_bed.hf$variable)

heatMapPlotter <- function (table, main="", weekNotation=F, legPos="right", mode="default", xlab="", ylab="")
{
  #Change weeks by Development and habituation notation
  if (weekNotation == T)
  {           
    table$week <- paste ("Dev Phase", table$period-1, sep = " ") 
    levels(table$week) <- c (levels(table$week), "Dev Phase")          
    angleY = 330
  }
  else
  {
    #only numbers on the y axis of the plot
    if (weekNotation == "N")
    {
      table$week <- table$period-1 
      levels(table$week) <- c (levels(table$week))
      angleY = 0
    }
    else
    {
      table$week <- paste ("week", table$period, sep = "_")
      angleY = 330
    }
  }
  
  #Checking mode for setting suitable color scale
  if (mode == "pvalues")
  {           
    colorsSc = c ('black', 'black', 'black', 'yellow', 'cyan', 'black','black',  'black')
    valuesSc   = c (-100,    -0.08,   -1.08,     -1,         0.00000000000000000001,         0.08,   0.08,    100)
    limitsSc= c (-0.06,0.06)
    breaksSc   = c (-0.05, -0.01, 0.01, 0.05)
    labelsSc = c (">0.05", "0.01", "0.01", ">0.05")
    legName = "p-value"          
  }
  else
  {
    colorsSc = c ('green', 'green', 'green', 'black', 'black', 'red', 'red', 'red')
    valuesSc   = c (-10,  -3, -3, 0, 0, 3, 3, 10)
    limitsSc= c (-3,3)
    breaksSc = c (-3, -2, -1, 0, 1, 2, 3)
    labelsSc = c ("<-3","-2","-1","0", "1", "2", ">3")
    legName = "Fold change"
  }
  
  table$period <- as.numeric (table$period)
  table$week <- with (table, reorder (week, period,))
  table$chVar <- paste (table$channel, table$varOrder, table$variable, sep = " ")

  (p <- ggplot (table, aes (week, chVar)) + geom_tile (aes (fill = foldChange),
                                                       colour = "white") + 
      geom_text(aes(label=stars), color="white", size=7) +
      scale_fill_gradientn (guide = "colorbar",
                            colours = colorsSc,
                            values = valuesSc,
                            limits = limitsSc,
                            breaks   = breaksSc,
                            labels = labelsSc,
                            name = legName,
                            rescaler = function(x,...) x,
                            oob = identity) + ggtitle(main) + theme (legend.position = "none"))#no legend
  base_size <- 11
  
  p + theme_grey (base_size = base_size) + labs (x = xlab,
                                                 y = ylab) + scale_x_discrete (expand = c(0, 0)) +
    scale_y_discrete (expand = c(0, 0), labels = table$variable) + 
    theme (axis.ticks = element_blank(),
           legend.position = legPos,
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.title.x =  element_text (size = base_size * 1.4, face = "bold"),
           axis.title.y =  element_text (size = base_size * 1.4, face = "bold", angle = 90),
           axis.text.x = element_text(size = base_size * 1.4, angle = angleY,
                                      hjust = 0, face = "bold", colour = "black"),                           
           legend.text = element_text (size = base_size * 1.2),      
           legend.title = element_text (size = base_size *1.2, face = "bold"),      
           plot.title = element_text (size=base_size * 1.5, face="bold"),
           axis.text.y = element_text (size = base_size * 1.4,hjust = 0, face = "bold", colour = "black"))                          
}

data.frame_bed.hf$stars <- ""
data.frame_bed.hf <- data.frame_bed.hf [with (data.frame_bed.hf, order (period, channel, orderOut) ),]
data.frame_bed.hf <- data.frame_bed.hf [data.frame_bed.hf$period < 10,]

name_file <- paste ("heatmap_behavior_week", ".", image_format, sep="")
plot_width <- 12
plot_height <- 10
setwd(pwd)
heatMapPlotter (data.frame_bed.hf, main="\n",  weekNotation="N", legPos="right", xlab="\nWeeks\n", ylab="\n")
ggsave (file=name_file, width=plot_width, height=plot_height, dpi=300)
