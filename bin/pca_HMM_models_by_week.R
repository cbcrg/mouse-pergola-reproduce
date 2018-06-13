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
### Mice data segmentation from chromHMM visualized using Gviz            ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

### Execution example
## Rscript plot_HMM_segmentation.R --path_bed_files="path_to_bed_files"
library("ggplot2")
library("reshape2")
library("FactoMineR")

{
  if("Gviz" %in% rownames(installed.packages(lib.loc="/gviz/time")) == TRUE) {
    library("Gviz", lib.loc="/gviz/time")
  }
  else {
    library("Gviz")
  }
}

## bed files to GRanges
# library("GenomicRanges")
# library("rtracklayer")

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
      density_heatmaps_tracking.R
      
      Arguments:
      --path_bed_files=path_bed_files - character
      --help                          - print this text
      
      Example:      
      ./plot_HMM_segmentation.R --path_bed_files=\"path_bed_files\" --image_format=\"image_format\" \n")
  
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
  if (is.null (argsL$path_bed_files))
  {
    stop ("[FATAL]: path to bed files not provided", stderr())
  }
  else
  {
    path_bed_files <- argsL$path_bed_files
  }
}

# Initial time to plot
{
  if (is.null (argsL$ini_time))
  {
    ini_time=0
    write("[WARNING]: Initial time to plot not set, default 0", stderr())
  }
  else
  {
    ini_time <- as.integer(argsL$ini_time)
  }
}

# End time to plot
{
  if (is.null (argsL$end_time))
  {
    end_time <- 0
    write("[WARNING]: Initial time to plot not set, default 0", stderr())
  }
  else
  {
    end_time <- as.integer(argsL$end_time)
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

#############################
## Read transition prob
path_to_trans_tbl <- "/Users/jespinosa/git/cbcrg/mouse-pergola-reproduce/results/chromHMM/transition_prob.tbl"
data.frame_prob_trans <- read.csv(path_to_trans_tbl, header=F, sep="\t")
# head (data.frame_prob_trans) #del
colnames(data.frame_prob_trans) <- c("week", "group", "id_mouse", "prob_mode", "State_ini", "State_end", "prob")

# data.frame_prob_trans_nice <- data.frame_prob_trans
# data.frame_prob_trans_nice$state_trans_week <- paste (data.frame_prob_trans$week, 
                                                      # data.frame_prob_trans_nice$State_ini, 
                                                      # data.frame_prob_trans_nice$State_end,
                                                      # sep="_")
# head(data.frame_prob_trans) #del
# rm(clean_df)
# clean_df <- subset (data.frame_prob_trans, 
                    # select = c(id_mouse, week, State_ini, State_end, prob))
# clean_df <- subset (clean_df , State_end != "empty")

# clean_df$id_mouse <- data.frame_prob_trans_nice$id_mouse
# clean_df$week <- data.frame_prob_trans_nice$week
# clean_df$State_ini <- data.frame_prob_trans_nice$State_ini
# clean_df$State_end <- data.frame_prob_trans_nice$State_end
# clean_df$prob <- data.frame_prob_trans_nice$prob

## dealing with NAs
## not anymore now simply empty
# data.frame_prob_trans[is.na(data.frame_prob_trans)] <- 10000
# data.frame_prob_trans <- data.frame_prob_trans[ !is.na(data.frame_prob_trans$prob), ]
# clean_df <- subset (data.frame_prob_trans_clean, prob != NA) # does not work
data.frame_prob_trans$tag <- paste(data.frame_prob_trans$State_ini, data.frame_prob_trans$State_end, sep="_")

data.frame_prob_trans_clean <- subset (data.frame_prob_trans, 
                                       select = c(id_mouse, week, tag, prob))

# clean_df <- data.frame_prob_trans_clean
# head(clean_df)

# length(clean_df [which(clean_df$id_mouse == 2), 1])

# clean_df [which(clean_df$State_end=="empty"),]$State_end <- 2
# head(data.frame_prob_trans_clean)
# length(data.frame_prob_trans_clean[,1])

# 3
# 7
# 10
# 12
# 17
# 18

# clean_df <- data.frame_prob_trans_clean
# clean_df <- subset (data.frame_prob_trans_clean, id_mouse ==3 & week ==1)
# clean_df
# clean_df <- clean_df [which(clean_df$id_mouse == 7 & clean_df$week == 0), ]

# dcast (clean_df, id_mouse ~ week + State_ini + State_end, value.var="prob")
                      
# example_dcast <- dcast (data.frame_prob_trans_clean, id_mouse ~ week + State_ini + State_end, value.var="prob")

#############################
## Read transition prob
path_to_em_tbl <- "/Users/jespinosa/git/cbcrg/mouse-pergola-reproduce/results/chromHMM/emission_prob.tbl"
data.frame_prob_em <- read.csv(path_to_em_tbl, header=F, sep="\t")
# head (data.frame_prob_em) #del
colnames(data.frame_prob_em) <- c("week", "group", "id_mouse", "prob_mode", "State", "Emission_n", "Emission", "binary", "prob")
data.frame_prob_em_filt <- data.frame_prob_em [which(data.frame_prob_em$binary==1), ]
data.frame_prob_em_filt$tag <- paste(data.frame_prob_em_filt$State, data.frame_prob_em_filt$Emission, sep="_bin_")
data.frame_prob_em_clean <- subset (data.frame_prob_em_filt, 
                                       select = c(id_mouse, week, tag, prob))

df.trans_and_em <- rbind (data.frame_prob_trans_clean, data.frame_prob_em_clean)

pca_df <- dcast (df.trans_and_em , id_mouse ~ week + tag, value.var="prob")
# install.packages("ca")
# library ("ca")
# 
# ca_res<- ca(pca_df)
# ca_res
# plot(ca_res)

# head(pca_df )
row.names (pca_df) <- pca_df$id_mouse
res <- PCA(pca_df[ , (2:length(pca_df[1,]))])

# Variance of PC1 and PC2
var_PC1 <- round (res$eig [1,2])
var_PC2 <- round (res$eig [2,2])

pca2plot <- as.data.frame (res$ind$coord)
pca2plot$id_mouse <- pca_df$id_mouse
pca_models <- ggplot(pca2plot, aes(x=-Dim.1, y=-Dim.2)) + geom_point() +
                   geom_text(aes(label=id_mouse),hjust=0, vjust=0)
pca_models

############
## BARPLOT
df.bars <- cbind (as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE)), names(res$var$coord[,1])[order(res$var$coord[,1]^2,decreasing=TRUE)])
df.bars_to_plot <- as.data.frame(df.bars)
df.bars_to_plot$index <- as.factor (df.bars_to_plot$V2)
class (df.bars_to_plot$V1)
df.bars_to_plot$value <- as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE))
df.bars_to_plot$index <- factor(df.bars_to_plot$index, levels = df.bars_to_plot$index[order(df.bars_to_plot$value, decreasing=TRUE)])

bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
  ylim (c(0, 80)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC1\n", x = "", y="Contribution in %\n") +
  #   ## Plots circle for slides thesis
  # scale_x_discrete (labels=c("Latency", "Gallagher", "% periphery", "Whishaw","Distance", "% NE", "Speed")) +
  #   ## Plots circle for slides thesis
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot

circle_plot <- as.data.frame (res$var$coord)
labels_v <- row.names(res$var$coord)
p_circle_plot <- ggplot(circle_plot) + geom_point()







circle_plot <- as.data.frame (res$var$coord)
labels_v <- row.names(res$var$coord)
which (circle_plot$Dim.1 < 0)

neg_labels <- labels_v [which (circle_plot$Dim.1 < 0)]
neg_positions <- circle_plot [which (circle_plot$Dim.1 < 0), c(1,2)]

# change positions for labels
# neg_positions [2,2] <- neg_positions [2,2] - 0.03 
# neg_positions [3,2] <- neg_positions [3,2] + 0
# neg_positions [4,2] <- neg_positions [4,2] - 0.02

pos_labels <- labels_v [which (circle_plot$Dim.1 >= 0)]
pos_positions <- circle_plot [which (circle_plot$Dim.1 >= 0), c(1,2)]

angle <- seq(-pi, pi, length = 50)
df.circle <- data.frame(x = sin(angle), y = cos(angle))

#aes(x=PC1, y=PC2, colour=gentreat )) 
p_circle_plot <- ggplot(circle_plot) + 
  # geom_segment (data=circle_plot, aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, color="red") +
  xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
  geom_text (data=neg_positions, aes (x=Dim.1, y=Dim.2, label=neg_labels, hjust=1.2), show.legend = FALSE, size=5) + 
  geom_text (data=pos_positions, aes (x=Dim.1, y=Dim.2, label=pos_labels, hjust=-0.3), show.legend = FALSE, size=5) +
  # geom_vline (xintercept = 0, linetype="dotted") +
  # geom_hline (yintercept=0, linetype="dotted") +
  labs (title = "PCA of the variables\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) #+
  # geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)
p_circle_plot
