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
### Creating bins from meal duration to train the HMM using ChromHMM      ###
#############################################################################

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
      distro_meals_to_bin.R
      
      Arguments:
      --path_bed_files=path_bed_files - character
      --n_bins=n_bins                 - integer
      --image_format=image_format     - character
      --help                          - print this text
      
      Example:
      ./distro_meals_to_bin.R --path_bed_files=\"path_bed_files\" --image_format=\"image_format\" \n")
  
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
  if (is.null (argsL$path_bed_files))
  {
    stop ("[FATAL]: Path to files is mandatory")
  }
  else
  {
    path_bed_files <- argsL$path_bed_files
  }
}

# path to files
{
  if (is.null (argsL$n_bins))
  {
    n_bins <- 3
    warning ("[Warning]: number of bins set to default (3)")
  }
  else
  {
    n_bins <- as.integer(argsL$n_bins)
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

## Loading libraries
library ("ggplot2")

## Loading parameters for plots
source("https://gist.githubusercontent.com/JoseEspinosa/307430167b05e408ac071b8724701abf/raw/06b26f764953ceea53d334190d7b736308e5141d/param_plot_publication.R")

# Whole experiment
# path_bed_files <- "/Users/jespinosa/git/mouse_chrom_hmm/results_bed"
pwd <- getwd()
setwd(path_bed_files)
files <- list.files(pattern=paste("tr_.*.bed$", sep=""))

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
  df$mouse <- mouse
  df$genotype <- genotype
  df$genotype <- factor(df$genotype, levels=c("wt", "KO"), labels=c("wt", "KO"))
  df$mouse <- mouse
  df$data_type <- data_type
  df$phase <- phase
  df$exp_phase <- gsub("_", " ", exp_phase)
  df$group2plot <- paste (phase, data_type)
  data.frame_bed <- rbind(data.frame_bed, df)
}

setwd(pwd)

data.frame_bed$length <- data.frame_bed$V3 - data.frame_bed$V2

## binning base on number of means
intercepts <- unname(quantile(data.frame_bed$length, probs = seq(1/n_bins, 1, 1/n_bins)))
intercepts <- intercepts[1:(length(intercepts)-1)]

write(cat(intercepts), stdout())

x <- data.frame_bed$length 
x_l <- "\nLength of meals (log seconds)"
y_l <- "Counts\n"
main_title <- "Distribution of meal length\n"

# log of length
data.frame_bed$length_log <- log(data.frame_bed$length)

ggplot(data.frame_bed, aes(x=length)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.1) +
  scale_x_continuous(breaks=sort(c(0, 10, 100, 1000, intercepts)), trans="log1p", expand=c(0,0)) +
  # scale_x_continuous(breaks=c(0, 10, 100, 1000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0, 25000, 5000), expand=c(0,0)) +
  geom_vline(xintercept=intercepts, linetype="dashed") +
  labs(title=main_title, x=x_l, y=y_l) +
  theme(plot.title = element_text(hjust = 0.5))# +

plot_width <- 14
plot_height <- 8

name_out <- paste ("meal_length_distro_binned", ".", image_format, sep="")
ggsave (file=name_out, width = plot_width, height=plot_height)
