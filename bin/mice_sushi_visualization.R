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
### Mice data set visualization using Sushi                               ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

library("ggplot2")
library("Sushi")
library("grid")

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
        --path_to_bedGraph_files=path_bedg_files - character
        --path_to_phases_file=path_phases_files  - character
        --help                                   - print this text
        
        Example:
        
        ./mice_sushi_visualization.R --f_experiment_info=\"path_to_file_experiment_info\" --path_bed_files=\"path_bed_files\" --path_to_bedGraph_files=\"path_bedg_files\" --path_to_phases_file=\"path_exp_phases\"  \n")
    
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
        write("[WARNING]: path_to_phases_file arg not provided", stderr())
        phases_file <- NULL
    }
    else
    {
        phases_file <- argsL$path_to_phases_file
    }
}

cb_palette <- c("#999999", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7", "#000000", 
                "#00009B")

color_heatmap <- 'blue'

#############################
## Read files bed files
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

## Assign a colo to each group
l_gr_color <- mapply(function(x, col) list(col),
                     levels(b2v$condition),
                     cb_palette[1:length(levels(b2v$condition))])

data_bed_events <- lapply(b2v$path, function (bed) { 
  
    name_id <- sub("tr_(.*?)_dt_.*\\.bed", "\\1", basename(bed))
    type_id <- sub("tr_.*_dt_(.*?)\\.bed", "\\1", basename(bed))
    
    bed_tbl <- read.csv(file=bed, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    bed_tbl$name <- paste(name_id, "_", type_id, sep="")
    bed_tbl$id <- name_id
    bed_tbl$color <- l_gr_color[[exp_info$condition [exp_info$sample == gsub("\\.bed", "", basename(bed))]]]
    
    return (bed_tbl)
})

#############################
## Read files bedGraph files
bedg_dir <- file.path(path_bedG_files)
perg_bedg_files <- sapply(exp_info$sample, function(id) file.path(bedg_dir, paste(id, ".bedGraph", sep="")))

bg2v<-b2v
bg2v <- dplyr::mutate(bg2v, path = perg_bedg_files)

max_value <- -1000000

data_bedg_win <- lapply(bg2v$path, function (bedg) { 
    
    name_id <- sub("tr_(.*?)_dt_.*\\.bedGraph", "\\1", basename(bedg))
    type_id <- sub("tr_.*_dt_(.*?)\\.bedGraph", "\\1", basename(bedg))    
    bedg_tbl <- read.csv(file=bedg, header=FALSE, sep="\t", stringsAsFactors=FALSE)
#     max_value <<- max (max_value, max(bedg_tbl$V4))
    bedg_tbl$name <- paste(name_id, "_", type_id, sep="")
    bedg_tbl$color <- l_gr_color[[exp_info$condition [exp_info$sample == gsub("\\.bedGraph", "", basename(bedg))]]]
    
    return (bedg_tbl)
})

## I set max_value to 0.5 to make the representation like in the Gviz case
max_value <- 0.5

## standardize values to a range of 0 to 1
unite_scale <- function (v, min=0, max=0.5) {    
    if (v > max) {return (1)}
    return ((v - min) / (max - min))
}

# parameters for sushi
chrom            = "chr1"
chromstart       = 0
chromend         = 100000
chromend         = 3635631

###########################################
## Read files bed phases files if available
{
    if (!is.null (phases_file)) {
        bed_phases <- read.csv(file=phases_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    }
    else {
        bed_phases <- data.frame(chrom=chrom, start=chromstart, end=chromend, V4=".", V5=1, V6="+", V7=chromstart, V8=chromend, V9)
    }
}

## Plot
pdf(paste("mice_sushi_viz", ".pdf", sep=""))

split.screen (c(2, 1)) 

## adding a n empty plots for title
n=3
# split.screen(c(length(data_bedg_win)+n, 1), screen = 1)
# split.screen(c(length(data_bed_events)+length(data_bedg_win)+n, 1), screen = 1)
split.screen(c(length(data_bed_events) + n, 1), screen = 1)
split.screen(c(length(data_bedg_win) + n, 1), screen = 2)

screen(1)

par(mar=c(0.1,1,2,0.1))

plot(1, type="n", axes=F, xlab="", ylab="")

# labelplot("A ", title, letteradj=-.025)

i=3+n
j=1

for (bed_i in seq_along(data_bed_events)) {
    screen( i )    
    par(mar=c(0.1, 2, 0.1, 2))
    plotBed(data_bed_events[[bed_i]], chrom, chromstart, chromend,
            row='supplied', color=data_bed_events[[bed_i]]$color[1])
    
    axis(side=2, lwd.tick=0, labels=FALSE, col=data_bed_events[[bed_i]]$color[1])
    mtext(data_bed_events[[bed_i]]$id[1], side=2, cex=0.5, col=data_bed_events[[bed_i]]$color[1])
    
    i=i+1
}

screen(2)

for (bedg_i in seq_along(data_bedg_win)) {
    screen(i)
    par(mar=c(0.1, 2, 0.1, 2))
    
    plotBed(data_bedg_win[[bedg_i]], chrom, chromstart, chromend, 
            row='supplied',  
            color=opaque(color_heatmap, transparency=unite_scale(data_bedg_win[[bedg_i]]$V4, max=max_value)))  
   
    axis(side=2, lwd.tick=0, labels=FALSE, col=data_bed_events[[bedg_i]]$color[1])

    mtext(data_bed_events[[bedg_i]]$id[1], side=2, cex=0.5, col=data_bed_events[[bedg_i]]$color[1])

    i=i+1
    j=j+1
}

screen(i)
par(mar=c(0.1, 2, 0.1, 2))
plotBed(bed_phases, chrom, chromstart, chromend, row='supplied')
mtext("Phases", side=2, cex=0.5, col="black")

min_heatmap <- 0
max_heatmap <- max_value
color_min <- 'white'
color_max <- color_heatmap

x <- runif(length(unique(exp_info$condition)),0,100)
y <- runif(length(unique(exp_info$condition)),100,200)

df_legend <- data.frame(x, y, gsub("_", " ", unique(exp_info$condition)))
colnames(df_legend) <- c("x", "y", "names")
df_legend$names <- ordered(gsub("_", " ", df_legend$names), levels = gsub("_", " ", levels(exp_info$condition)))
color_by_tr <- unlist(l_gr_color[levels(exp_info$condition)])
names(color_by_tr) <- gsub("_", " ", levels(exp_info$condition))
size_text_leg <- 8
df_empty <- data.frame()

plot_legends <- ggplot(df_empty) + geom_point() + theme(panel.border = element_blank(), panel.background = element_blank())
size_box_leg <- 3
plot_legends <- plot_legends + geom_point(data=df_legend, aes(x=x, y=y, colour = names), shape=15, size=size_box_leg) +
    scale_colour_manual (values=color_by_tr) + 
    guides(color=guide_legend(title=NULL)) + 
    theme(legend.position="bottom", legend.justification=c(0, 0), 
          legend.text=element_text(size=size_text_leg),
          legend.key = element_rect(fill = "white", colour = "white")) + 
    geom_blank()

## Adding heatmap scale to the legend
bedGraphRange <- c(0, 0.5)
plot_legends <- plot_legends + geom_point(data=df_legend, aes(x=x, y=y, fill = 0)) +
    scale_fill_gradientn (guide = "colorbar",
                          colours = c(color_min, color_max),
                          values = c(bedGraphRange[1], bedGraphRange[2]),
                          limits = c(bedGraphRange[1], bedGraphRange[2]),
                          breaks = c(bedGraphRange[1], bedGraphRange[2]),
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

grid.draw (leg_groups)

dev.off()
