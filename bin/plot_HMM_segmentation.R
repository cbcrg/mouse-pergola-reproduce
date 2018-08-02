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

{
  if("Gviz" %in% rownames(installed.packages(lib.loc="/gviz/time")) == TRUE) {
    library("Gviz", lib.loc="/gviz/time")
  }
  else {
    library("Gviz")
  }
}

## bed files to GRanges
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
      --path_bed_files=path_bed_files - character
      --ini_time=ini_time             - integer
      --end_time=end_time             - integer
      --path_to_phases_file=path      - character
      --path_to_exp_phases_file=path  - character
      --image_format=image_format     - character
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
    write("[WARNING]: End time to plot not set, default 0", stderr())
  }
  else
  {
    end_time <- as.integer(argsL$end_time)
  }
}

# phases file
{
    if (is.null (argsL$path_to_phases_file))
    {
        phases_file <- ""
        write ("[WARNING]: path_to_phases_file arg not provided")
    }
    else
    {
        phases_file <- argsL$path_to_phases_file
    }
}

# experimental phases file
{
  if (is.null (argsL$path_to_exp_phases_file))
  {
      exp_phases_file <- ""
      write ("[WARNING]: path_to_exp_phases_file arg not provided")
  }
  else
  {
    exp_phases_file <- argsL$path_to_exp_phases_file
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
## Read files bed files
# end_time<-0
# ini_time <-0
bed_files <- list.files(path=path_bed_files, pattern="dense\\.bed$", full.names=TRUE)

bed_tracks <- lapply(bed_files, function (bed) {
  id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bed)
  bed_GR <- import(bed, format = "bed")
  tr <- AnnotationTrack(bed_GR, name = paste ("", id, sep=""),
                        fill=bed_GR$itemRgb)
  
  names(tr) <- id
  
  return (tr)
} )

bed_GR <- import (bed_files[[1]], format = "bed")
states_names <- unique (bed_GR$name)
color_by_tr <- unique (bed_GR$itemRgb)
n_states <- length (states_names)
names (color_by_tr) <- states_names
state_names <- c("Long meals ", "Regular meals (2/3 min) ", "Inactive ", "Short meals ")


## Reorder to show plots in the same order
phases_file <- file.path(phases_file)
name_phases_tr <- ""
phases_color <- 'gray'
col_back_title="brown"

{
  if(file.exists(phases_file)) {
    bed_phases <- import(phases_file, format = "BED")
    phases_tr <- AnnotationTrack(bed_phases, name = paste ("", name_phases_tr, sep=""),
                                 fill = c("#EADAFF"),
                                 background.title = col_back_title, col=NULL)
  }
  else {
    phases_tr <- NULL
  }
}

exp_phases_file <- file.path(exp_phases_file)
name_exp_phases_tr <- ""
phases_color <- 'gray'
col_back_title="brown"

{
  if(file.exists(exp_phases_file)) {
    bed_exp_phases <- import(exp_phases_file, format = "BED")
    exp_phases_tr <- AnnotationTrack(bed_exp_phases, name = paste ("", name_exp_phases_tr, sep=""),
                                 fill = c("#999999", "#000000"),
                                 background.title = col_back_title, col=NULL)
  }
  else {
    exp_phases_tr <- NULL
  }
}




names(bed_tracks) <- as.numeric(gsub(".+tr_(\\d+)(_.+$)", "\\1", bed_files))
id_mice <- sort(as.numeric(gsub(".+tr_(\\d+)(_.+$)", "\\1", bed_files)))

mice_id_odd <- id_mice [id_mice %% 2 != 0]
mice_id_even <- id_mice [id_mice %% 2 == 0]

bed_tracks <- bed_tracks[as.character(c(mice_id_odd, mice_id_even))]

## Plot

size_labels <- 16
cex_gtrack <- 1.4
g_tr <- GenomeAxisTrack()

## Empty tracks for placing legend
ctracks <- list()
for (i in 1:2) {
  ctracks[[i]] <- CustomTrack(plottingFunction=function(GdObject, prepare, ...) {
    if(!prepare) grid.text("")
    return(invisible(GdObject))
  }, variables=list(i=i))
  displayPars(ctracks[[i]]) <- list(background.title="transparent")
}

## creating a legend
x <- runif(c(1: n_states), 0, 100)
y <- runif(c(1: n_states), 100, 200)

df_legend <- data.frame(x, y, as.character(c(1: n_states)))
colnames(df_legend) <- c("x", "y", "n_states")

size_text_leg <- 18
# size_text_leg <- 10
df_empty <- data.frame()

plot_legends <- ggplot(df_empty) + geom_point() + 
                theme(panel.border = element_blank(),
                panel.background = element_blank())
size_box_leg <- 6
# size_box_leg <- 4

{
  if (image_format == 'pdf') {
      size_box_leg <- size_box_leg * 3
      size_text_leg <- size_text_leg * 3
      cex_gtrack <- cex_gtrack * 2
      size_labels <- size_labels * 2
  }
}

plot_legends <- plot_legends + geom_point(data=df_legend, 
                                          aes(x=x, y=y, colour = n_states), 
                                          shape=15, size=size_box_leg) +
                scale_colour_manual (values=color_by_tr, labels=state_names) +
                guides(color=guide_legend(title=NULL)) +
                theme(legend.position="bottom", legend.justification=c(1, 0),
                      legend.text=element_text(size=size_text_leg),
                      legend.key = element_rect(fill = "white", colour = "white")) +
                geom_blank()

## Extract Legend 
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

leg_groups <- g_legend(plot_legends) 

plot_name <- "segmentation_HMM"

{
  if (image_format == 'tiff' | image_format == 'tif') {
    tiff(paste(plot_name, ".", image_format, sep=""), width = 80 , height = 40, units = "cm", res=300)
  }
  else if (image_format == 'pdf') {
    pdf(paste(plot_name, ".", image_format, sep=""), height=40, width=80)
  }
  else if (image_format == 'png') {        
    png(paste(plot_name, ".", image_format, sep=""),  width = 80 , height = 40, units = "cm", res=300)
  }
  else {
    stop (paste("Unknow image file format:", image_format, sep=" "))
  }
}

cex_gtrack <- 1
v_height_tracks <- c(0.2, rep(0.1, length(unlist(bed_tracks))), 0.1, 0.1, rep(0.1,length(unlist(ctracks))))

{ 
  if (end_time == 0) {
    plotTracks(c(g_tr, unlist(bed_tracks), unlist(phases_tr), unlist(exp_phases_tr), unlist(ctracks)),
    # plotTracks(c(g_tr, unlist(bed_tracks)),
               from=ini_time, 
               stacking="dense", 
               collapse=FALSE, 
               shape = "box", 
               col=NULL,
               fontsize=size_labels, 
               cex=cex_gtrack,
               sizes=v_height_tracks)
  }
  else {
    plotTracks(c(g_tr, unlist(bed_tracks), unlist(phases_tr), unlist(exp_phases_tr), unlist(ctracks)),
    # plotTracks(c(g_tr, unlist(bed_tracks)),
               from=ini_time, to=end_time,
               stacking="dense", 
               collapse=FALSE, 
               shape = "box", 
               col=NULL,
               stacking = "dense",
               fontsize=size_labels, 
               cex=cex_gtrack,
               sizes=v_height_tracks)
  }
}

grid.draw (leg_groups)

dev.off()
