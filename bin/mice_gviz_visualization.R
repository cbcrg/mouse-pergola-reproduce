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
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. Aug 2017                   ###
#############################################################################
### Mice data set visualization using Gviz                                ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

### Execution example
## Rscript density_heatmaps_tracking.R --files_str1="list_files_str1" --files_str2="list_files_str2"
library(ggplot2)
# library("Biostrings")
# library("devtools")
library("Gviz")

## bedfiles to GRanges
library("GenomicRanges")
library("rtracklayer")

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
      --f_experiment_info=experiment_info      - character
      --path_bed_files=path_to_bed_files       - character
      --path_bedg_files=path_to_bedGraph_files - character
      --help                           - print this text

      Example:
      ./plot_speed_motion_mean.R --f_experiment_info=\"path_to_file_experiment_info\" --path_bed_files=\"path_bed_files\" --path_to_bedGraph_files=\"path_bedg_files\" \n")

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
  if (is.null (argsL$f_experiment_info))
  {
    stop ("[FATAL]: f_experiment_info arg is mandatory")
  }
  else
  {
    experiment_info <- argsL$f_experiment_info
  }
}

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

{
  if (is.null (argsL$path_to_bedGraph_files))
  {
    stop ("[FATAL]: path to bedGraph files not provided", stderr())
  }
  else
  {
    path_bedG_files <- argsL$path_to_bedGraph_files
  }
}

{
  if (is.null (argsL$path_to_phases_file))
  {
    stop ("[WARNING]: path_to_phases_file arg not provided")
  }
  else
  {
    phases_file <- argsL$path_to_phases_file
  }
}

#############################
## Read files bedGraph files
b2v <- exp_info <- read.table(file.path(experiment_info), header = TRUE, stringsAsFactors=FALSE)
exp_info$condition <- as.factor(exp_info$condition)
exp_info$condition <- ordered(exp_info$condition, levels = c("WT_food_sc",
                                                             "WT_nic_food_sc", 
                                                             "CB1_food_sc", 
                                                             "CB1_nic_food_sc",
                                                             "WT_food_fat",
                                                             "WT_nic_food_fat",
                                                             "CB1_food_fat",                                                             
                                                             "CB1_nic_food_fat"))
b2v <- exp_info

bed_dir <- file.path(path_bed_files)

{
  if (length(exp_info$sample) != length(unique(exp_info$sample))) {
    stop ("Sample names duplicated in configuration file")}
}

perg_bed_files <- sapply(exp_info$sample, function(id) file.path(bed_dir, paste(id, ".bed", sep="")))

b2v <- dplyr::mutate(b2v, path = perg_bed_files, header = TRUE, stringsAsFactors=FALSE)
bedg_dir <- file.path(path_bedG_files)
perg_bedg_files <- sapply(exp_info$sample, function(id) file.path(bedg_dir, paste(id, ".bedGraph", sep="")))

bg2v<-b2v
bg2v <- dplyr::mutate(bg2v, path = perg_bedg_files)

cb_palette <- c("#999999", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7", "#000000", 
                "#00009B")

l_gr_color <- mapply(function(x, col) list(col),
                     levels(b2v$condition),
                     cb_palette[1:length(levels(b2v$condition))])

bed2pergViz <- function (data_df, gr_df, format_f="BED") {
  grps <- as.character(gr_df[[setdiff(colnames(gr_df), 'sample')]])
  r <- lapply(unique(grps),
              function(g) {
                gr_samps <- grps %in% g
                gr_files <- data_df$path[gr_samps]
                lapply(gr_files, function (bed) {
                  id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bed)
                  bed_GR <- import(bed, format = format_f)

                  if (format_f == "BED") {
                    tr <- AnnotationTrack(bed_GR, name = paste ("", id, sep=""),
                                          fill=l_gr_color[[g]],
                                          background.title = l_gr_color[[g]], col=NULL)
                  }

                  if (format_f == "bedGraph") {
                    min_v <<- floor (min(bed_GR$score))
                    max_v <<- ceiling (max(bed_GR$score))
                    tr <- bed_GR
                  }

                  return (tr) })
              })

  names(r) <- unique(grps)
  return (r)
}

l_gr_annotation_tr_bed <- bed2pergViz (b2v, exp_info)

## Reorder to show plots in the same order
phases_file <- file.path(phases_file)
name_phases_tr <- "Phases"
phases_color <- 'gray'
col_back_title="brown"

{ 
  if(file.exists(phases_file)) {
    bed_phases <- import(phases_file, format = "BED")
    phases_tr <- AnnotationTrack(bed_phases, name = paste ("", name_phases_tr, sep=""),
                                 fill = phases_color,
                                 background.title = col_back_title, col=NULL)    
  }
  else {
    phases_tr <- NULL    
  }
}

l_granges_bg <- bed2pergViz (bg2v, exp_info, "bedGraph") 

min_heatmap <- 0
max_heatmap <- 0.5
color_min <- 'white'
color_max <- 'blue'

# to make the order equal in both panels it should be reorder before applying colors (numeric ordering)
l_granges_bg <- l_granges_bg[levels(exp_info$condition)]

l_gr_data_tr_bg <- lapply (seq_along(l_granges_bg), function (i_group_exp) {
  lapply (seq_along (l_granges_bg[[i_group_exp]]),  function (i_track) {                                                     
    granges_obj <-l_granges_bg[[i_group_exp]][[i_track]]
    tr_name <- names(l_granges_bg[[i_group_exp]][i_track])                                                     
    id <- gsub("^tr_(\\d+)(_dt.*$)", "\\1", tr_name)                                                     
    d_track <- DataTrack(granges_obj,
                         type="heatmap", ylim = c(min_heatmap, max_heatmap),
                         background.title = l_gr_color[[i_group_exp]],
                         gradient=c(color_min, color_max), 
                         showAxis = F, name = id)
    return (d_track)
  })
})

names (l_gr_data_tr_bg) <- names(l_granges_bg)

## reorder groups
l_gr_annotation_tr_bed <- l_gr_annotation_tr_bed[levels(exp_info$condition)]
l_gr_data_tr_bg <- l_gr_data_tr_bg [levels(exp_info$condition)]

size_labels <- 12
cex_gtrack <- 0.7
g_tr <- GenomeAxisTrack()

## creating a legend
x <- runif(length(unique(exp_info$condition)),0,100)
y <- runif(length(unique(exp_info$condition)),100,200)
# df_legend <- data.frame(x, y, ordered(df_legend$names, levels = levels(exp_info$condition)))
df_legend <- data.frame(x, y, gsub("_", " ", unique(exp_info$condition)))
colnames(df_legend) <- c("x", "y", "names")
df_legend$names <- ordered(gsub("_", " ", df_legend$names), levels = gsub("_", " ", levels(exp_info$condition)))
color_by_tr <- unlist(l_gr_color[levels(exp_info$condition)])
names(color_by_tr) <- gsub("_", " ", levels(exp_info$condition))
size_text_leg <- 14
df_empty <- data.frame()

plot_legends <- ggplot(df_empty) + geom_point() + 
                theme(panel.border = element_blank(), panel.background = element_blank())
size_box_leg <- 5
plot_legends <- plot_legends + geom_point(data=df_legend, aes(x=x, y=y, colour = names), shape=15, size=size_box_leg) +
                scale_colour_manual (values=color_by_tr) + 
                guides(color=guide_legend(title=NULL)) + 
                theme(legend.position="bottom", legend.justification=c(0, 0), 
                      legend.text=element_text(size=size_text_leg),
                      legend.key = element_rect(fill = "white", colour = "white")) + 
                geom_blank()

## Adding heatmap scale to the legend
bedGraphRange <- c(0,0.5)
plot_legends <- plot_legends + geom_point(data=df_legend, aes(x=x, y=y, fill = 0)) +
  scale_fill_gradientn (guide = "colorbar",
                        colours = c(color_min, color_max),
                        values = c(bedGraphRange[1], bedGraphRange[2]),
                        limits = c(bedGraphRange[1], bedGraphRange[2]),
                        breaks   = c(bedGraphRange[1], bedGraphRange[2]),
                        labels = c(bedGraphRange[1], paste(bedGraphRange[2],"    ", sep="")),
                        name = "",
                        rescaler = function(x,...) x,                                        
                        oob = identity) + theme (legend.position = "none") + 
  theme(legend.position="bottom", legend.justification=c(1,0), legend.text=element_text(size=size_text_leg)) +
  geom_blank()  

## Extract Legend 
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

leg_groups <- g_legend(plot_legends) 

## Empty tracks for placing legend
ctracks <- list()
for (i in 1:10) {
  ctracks[[i]] <- CustomTrack(plottingFunction=function(GdObject, prepare, ...) {
    if(!prepare) grid.text("")
    return(invisible(GdObject))
  }, variables=list(i=i))
  displayPars(ctracks[[i]]) <- list(background.title="transparent")
}

name_file <- "mice_gviz_viz.tiff"
# name_file <- "mice_gviz_viz.png"
# png(name_file, width = 45 , height = 34, units = "cm", res=300)
# png(name_file, width = 2000 , height = 1800, res=100)
tiff(name_file, width = 45 , height = 34, units = "cm", res=300)
# names(l_gr_annotation_tr_bed)
p <- plotTracks(c(g_tr, unlist(l_gr_annotation_tr_bed), unlist(l_gr_data_tr_bg),  phases_tr, unlist(ctracks)),
           from=0, to=3628800, 
           ylim=c(0,0.5),                                                     
           shape = "box", stacking = "dense",
           fontsize=size_labels, cex=cex_gtrack)
grid.draw (leg_groups)
dev.off()
