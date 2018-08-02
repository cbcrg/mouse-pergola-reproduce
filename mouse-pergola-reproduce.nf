#!/usr/bin/env nextflow

/*
 *  Copyright (c) 2014-2018, Centre for Genomic Regulation (CRG).
 *  Copyright (c) 2014-2018, Jose Espinosa-Carrasco and the respective authors.
 *
 *  This file is part of Pergola.
 *
 *  Pergola is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Pergola is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Jose Espinosa-Carrasco. CB-CRG. March 2017
 *
 * Script to reproduce Pergola paper figures of CB1 mice experiment
 */

params.recordings     = "$baseDir/small_data/mouse_recordings/"
params.mappings       = "$baseDir/small_data/mappings/b2p.txt"
params.mappings_bed   = "$baseDir/small_data/mappings/bed2pergola.txt"
params.phases         = "$baseDir/small_data/phases/exp_phases.csv"
params.mappings_phase = "$baseDir/small_data/mappings/f2g.txt"
params.exp_info       = "$baseDir/small_data/mappings/exp_info.txt"
params.tbl_chromHMM   = "$baseDir/small_data/chromHMM_files/cellmarkfiletable"
params.n_bins_HMM     = 5
params.n_states_HMM   = 4
params.output         = "files/"
params.image_format   = "png"

log.info "Mouse - Pergola - Reproduce  -  version 0.2.0"
log.info "====================================="
log.info "mice recordings          : ${params.recordings}"
log.info "mappings                 : ${params.mappings}"
log.info "mappings bed             : ${params.mappings_bed}"
log.info "experimental phases      : ${params.phases}"
log.info "mappings phases          : ${params.mappings_phase}"
log.info "experimental info        : ${params.exp_info}"
log.info "chromHMM config table    : ${params.tbl_chromHMM}"
log.info "HMM number of bins       : ${params.n_bins_HMM}"
log.info "HMM number of states     : ${params.n_states_HMM}"
log.info "output                   : ${params.output}"
log.info "image format             : ${params.image_format}"
log.info "\n"

// Example command to run the script
/*
nextflow run mouse-pergola-reproduce.nf \
  --recordings='small_data/mouse_recordings/' \
  --mappings='small_data/mappings/b2p.txt' \
  --mappings_bed='small_data/mappings/bed2pergola.txt' \
  --phases='small_data/phases/exp_phases.csv' \
  --mappings_phase='small_data/mappings/f2g.txt' \
  --exp_info='small_data/mappings/exp_info.txt' \
  --tbl_chromHMM="small_data/chromHMM_files/cellmarkfiletable" \
  --n_bins_HMM=5 \
  --n_states_HMM=4 \
  --image_format='png' \
  -with-docker
*/

/*
 * Input parameters validation
 */
mapping_file = file(params.mappings)
mapping_bed_file = file(params.mappings_bed)
mapping_file_bG = file(params.mappings)
mapping_file_phase = file(params.mappings_phase)

exp_phases = file(params.phases)
exp_info = file(params.exp_info)

/*
 * Input files validation
 */
if( !mapping_file.exists() ) exit 1, "Missing mapping file: ${mapping_file}"
if( !mapping_file_phase.exists() ) exit 1, "Missing mapping phases file: ${mapping_file_phase}"
if( !exp_phases.exists() ) exit 1, "Missing phases file: ${exp_phases}"
if( !exp_info.exists() ) exit 1, "Missing experimental info file: ${exp_info}"
cell_mark_file_tbl = Channel.fromPath(params.tbl_chromHMM)

// HMM parametrization
n_bins = params.n_bins_HMM
n_states = params.n_states_HMM

/*
 * Read image format
 */
image_format = "${params.image_format}"

if( image_format.matches('tiff') ) {
    println "WARNING: Deeptools figures will be created in png format as tiff format is not available.\n"
    image_format_deeptools = 'png'
}
else {
    image_format_deeptools = image_format
}

/*
 * Create a channel for mice recordings
 */
Channel
    .fromPath( "${params.recordings}intake*.csv" )
    .ifEmpty { error "Cannot find any CSV file with mice data" }
    .set { mice_files }

mice_files.into { mice_files_bed; mice_files_bedGraph }

/*
 * Create a channel for mice recordings
 */
Channel
    .fromPath( params.recordings )
    .set { mice_files_preference }

/*
 * Calculates mice preference statistics
 */
process behavior_by_week {

  	input:
  	file file_preferences from mice_files_preference
  	file mapping_file
    file mapping_file_phase
    file exp_phases

  	output:
  	file 'behaviors_by_week' into d_behaviors_by_week
  	stdout into max_time
    file 'exp_phases' into exp_phases_bed_to_wr, exp_phases_bed_to_wr2, exp_phases_bed_to_fraction, exp_phases_bed_to_hmm

	file 'stats_by_phase/phases_dark.bed' into exp_circadian_phases_sushi, exp_circadian_phases_gviz, days_bed_igv, days_bed_shiny, days_bed_deepTools, days_bed_gviz_hmm
    file 'Habituation_dark.bed' into bed_dark_habituation, bed_dark_habituation_groups
    file 'Development_dark.bed' into bed_dark_development, bed_dark_development_groups
    file 'whole_experiment_dark.bed' into whole_experiment_dark
    file 'whole_experiment_light.bed' into whole_experiment_light

  	"""
  	mice_stats_by_week.py -f "${file_preferences}"/intake*.csv -m ${mapping_file} -s "sum" -b feeding -p ${exp_phases} \
  	                      -mp ${mapping_file_phase}

  	mkdir stats_by_phase
  	mkdir behaviors_by_week
  	cp exp_phases.bed exp_phases
  	tail +2 exp_phases > exp_phases_sushi
  	mv *.bed  stats_by_phase/
  	mv stats_by_phase/Development_dark.bed ./
  	mv stats_by_phase/Habituation_dark.bed ./
  	mv stats_by_phase/whole_experiment_dark.bed ./
  	mv stats_by_phase/whole_experiment_light.bed ./
  	mv *.tbl  behaviors_by_week/
  	"""
}

/*
 * Creates a heatmap that compares mice feeding behavior of high-fat mice against their controls
 */
process heatmap {

    publishDir "${params.output_res}/heatmap/", mode: 'copy', overwrite: 'true'

    input:
    file behaviors_by_week from d_behaviors_by_week

    output:
    file "*.${image_format}" into heatmap_behaviors

    """
    heatmap_behavior.R --path2files=${behaviors_by_week} --image_format=${image_format}
    """
}

def igv_files_by_group ( file ) {

    def map_id_group = [ "ctrl" : [1,3,5,7,9,11,13,15,17],
                         "hf" : [2,4,8,10,12,14,16,18] ]

    def id = file.split("\\_")[1]

    def food = file.split("\\_")[4].split("\\.")[0]

    def ext = file.split("\\.")[1]

    if ( map_id_group.get("ctrl").contains(id.toInteger()) && food == "sc" )
    return "igv/1_ctrl_food_sc/${id}.${ext}"

    if ( map_id_group.get("hf").contains(id.toInteger()) && food == "fat" )
    return  "igv/2_hf_food_fat/${id}.${ext}"

}

/*
 * Converts input behavioral trajectory of mice into BED files (discrete intervals)
 */
process convert_bed {

    publishDir params.output_res, mode: 'copy', pattern: "tr*food*.bed", saveAs: this.&igv_files_by_group

  	input:
  	file ('batch') from mice_files_bed
  	file mapping_file
  	file mapping_bed_file

  	output:
  	file 'tr*food*.bed' into bed_out, bed_out_shiny_p, bed_out_gviz, bed_out_sushi
  	file 'dir_bed' into dir_bed_to_bin
  	file 'phases_light.bed' into phases_night
  	file '*.fa' into out_fasta

  	"""
  	pergola_rules.py -i ${batch} -m ${mapping_file} -f bed -nt -e -bl -dl food_sc food_fat -d all

    shopt -s nullglob

	## ctrl
  	for f in {1,3,5,7,9,11,13,15,17}
  	do
  	    mkdir -p work_dir
  	    mkdir -p dir_bed_ctrl

  	    files=( tr_"\$f"_* )

  	    if (( \${#files[@]} )); then
  	        cd work_dir
  	        track_int=`ls ../"tr_"\$f"_"*`
  	        mv \${track_int} track_int
  	        echo -e "food_sc\tblack" > dict_color
  	        echo -e "food_fat\tblack" >> dict_color
  	        pergola_rules.py -i track_int -m ../${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color' -dl food_sc food_fat -d all
  	        in_f_sc=`ls tr_chr1*food_sc.bed`
            mv "\$in_f_sc" "`echo \$in_f_sc | sed s/chr1/\${f}/`"
  	        cd ..
  	        cp work_dir/tr*.bed ./dir_bed_ctrl/
  	        mv work_dir/tr*.bed ./
  	        mv work_dir/*.fa ./
        fi
  	done

    for f in {2,4,8,10,12,14,16,18}
  	do
  	    mkdir -p work_dir
        mkdir -p dir_bed_hf

  	    files=( tr_"\$f"_* )

  	    if (( \${#files[@]} )); then
  	        cd work_dir
  	        track_int=`ls ../"tr_"\$f"_"*`
  	        mv \${track_int} track_int
  	        echo -e "food_sc\torange" > dict_color
  	        echo -e "food_fat\torange" >> dict_color
  	        pergola_rules.py -i track_int -m ../${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color' -dl food_sc food_fat -d all
  	        in_f_sc=`ls tr_chr1*food_sc.bed`
            mv "\$in_f_sc" "`echo \$in_f_sc | sed s/chr1/\${f}/`"
  	        cd ..
  	        cp work_dir/tr*.bed ./dir_bed_hf/
  	        mv work_dir/tr*.bed ./
  	        mv work_dir/*.fa ./
        fi
  	done

    mkdir dir_bed

  	cp tr*.bed ./dir_bed
  	"""
}

result_dir_shiny_p = file("$baseDir/files")

result_dir_shiny_p.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_shiny_p"
}

bed_out_shiny_p.flatten().subscribe {
    it.copyTo( result_dir_shiny_p.resolve (  ) )
}

result_dir_IGV = file("results/igv/")


longest_fasta = out_fasta
                   .max { it.size() }

longest_fasta.subscribe {
    fasta_file = it
    fasta_file.copyTo ( result_dir_IGV.resolve ( "mice.fa" ) )
}

days_bed_shiny.subscribe {
    phases_file = it
    phases_file.copyTo ( result_dir_shiny_p.resolve ( "phases_dark.bed" ) )
}

/*
 * Converts input behavioral trajectory of mice into bedGraph files showing a continuous score along time windows (30 min)
 */
process convert_bedGraph {

    publishDir params.output_res, mode: 'copy', pattern: "tr*food*.bedGraph", saveAs: this.&igv_files_by_group

  	input:
    file ('batch_bg') from mice_files_bedGraph
  	file mapping_file_bG
  	val max from max_time.first()

  	output:
  	file '{tr_[1-9]_dt_*food*.bedGraph,tr_[1-9][0-9]_dt_*food*.bedGraph}' into bedGraph_out, bedGraph_out_shiny_p, bedGraph_out_gviz, bedGraph_out_sushi, bedGraph_out_bigwig
  	file '{tr_[1][1,3,5,7]_dt_*food*.bedGraph,tr_[1,3,5,7,9]_dt_*food*.bedGraph}' into bedGraph_out_bigwig_ctrl
  	file '{tr_[1][0,2,4,6,8]_dt_*food*.bedGraph,tr_[2,4,8]_dt_*food*.bedGraph}' into bedGraph_out_bigwig_hf
  	file 'chrom.sizes' into chrom_sizes,chrom_sizes_ctrl,chrom_sizes_hf, chrom_sizes_gr, chrom_sizes_chromHMM_binarize, chrom_sizes_chromHMM_l

    //file 'tr_10_12_14_16_18_2_4_8_dt_food_fat_food_sc.bedGraph' into bedGraph_hf
    //file 'tr_11_13_15_17_1_3_5_7_9_dt_food_sc.bedGraph' into bedGraph_ctrl

  	"""
  	pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt -e -dl food_sc food_fat -d all

  	#pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt \
  	#                 -e -dl food_sc food_fat -d all  -a join_all -t 1 3 5 7 9 11 13 15 17

  	#pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt \
  	#                 -e -dl food_sc food_fat -d all  -a join_all -t 2 4 8 10 12 14 16 18

  	# awk '{printf "%i", \$2/3600/24}' chrom.sizes
  	awk '{printf "%i", \$2/604800}' chrom.sizes
  	"""
}

bedGraph_out_shiny_p.flatten().subscribe {
    it.copyTo( result_dir_shiny_p.resolve ( ) )
}

exp_phases_bed_to_wr.subscribe {
    it.copyTo( result_dir_IGV.resolve ( 'exp_phases.bed' ) )
}

days_bed_igv.subscribe {
    it.copyTo( result_dir_IGV.resolve ( 'days_nights.bed' ) )
}

exp_phases_bed_to_wr2.subscribe {
    it.copyTo( result_dir_shiny_p.resolve ( 'exp_phases.bed' ) )
}

/*
 * Renders BED and bedGraph files using Gviz
 */
process gviz_visualization {

    publishDir "${params.output_res}/gviz", mode: 'copy', overwrite: 'true'

    input:
    file 'exp_info' from exp_info
    file 'bed_dir/*' from bed_out_gviz.collect()
    file 'bedgr_dir/*' from bedGraph_out_gviz.collect()
    file exp_phases_bed from exp_circadian_phases_gviz

    output:
    file "*.${image_format}" into gviz

  	"""
    mice_gviz_visualization.R --f_experiment_info=${exp_info} \
        --path_bed_files=bed_dir \
        --path_to_bedGraph_files=bedgr_dir \
        --path_to_phases_file=${exp_phases_bed} \
        --image_format=${image_format}
  	"""
}

/*
 * Renders BED and bedGraph files using Sushi
 */
process sushi_visualization {

    publishDir "${params.output_res}/sushi", mode: 'copy', overwrite: 'true'

    input:
    file 'exp_info' from exp_info
    file 'bed_dir/*' from bed_out_sushi.collect()
    file 'bedgr_dir/*' from bedGraph_out_sushi.collect()
    file exp_phases_bed from exp_circadian_phases_sushi

    output:
    file "*.${image_format}" into sushi

  	"""
    mice_sushi_visualization.R --f_experiment_info=${exp_info} \
        --path_bed_files=bed_dir \
        --path_to_bedGraph_files=bedgr_dir \
        --path_to_phases_file=${exp_phases_bed} \
        --image_format=${image_format}
  	"""
}

bedGraph_out_bigwig_short_name = bedGraph_out_bigwig.flatten().map {
    def content = it
    def name = it.baseName.replaceAll('tr_','').replaceAll('_dt_food_fat_food_sc','').replaceAll('_dt_food_sc','')
    [ content, name ]
}


/*
 * Convert bedGraph to bigWig files (deeptools input data)
 */
process bedgraph_to_bigWig {

  	input:
  	set file (bedgr_file), val (name) from bedGraph_out_bigwig_short_name
    file chrom_sizes from chrom_sizes.first()

    output:
    file '*.bw' into bigWig_matrix

    """
  	head -n -2 ${bedgr_file} > ${name}.trimmed

    bedGraphToBigWig ${name}.trimmed ${chrom_sizes} ${name}".bw"
    """
}

process bedgraph_ctrl_to_bigWig {

  	input:
  	file (bedgr_file) from bedGraph_out_bigwig_ctrl.flatten()
    file chrom_sizes from chrom_sizes_ctrl.first()

    output:
    file '*.bw' into bigWig_ctrl

    """
    bedGraphToBigWig ${bedgr_file} ${chrom_sizes} ${bedgr_file}".bw"
    """
}

process bedgraph_ctrl_to_bigWig {

  	input:
  	file (bedgr_file) from bedGraph_out_bigwig_hf.flatten()
    file chrom_sizes from chrom_sizes_hf.first()

    output:
    file '*.bw' into bigWig_hf

    """
    bedGraphToBigWig ${bedgr_file} ${chrom_sizes} ${bedgr_file}".bw"
    """
}

process bedgraph_to_mean_gr_bigWig {

  	input:
  	file bigwig_ctrl from bigWig_ctrl.collect()
    file bigwig_hf from bigWig_hf.collect()

    file chrom_sizes from chrom_sizes_gr.first()

    output:
    file 'ctrl.bw' into bigWig_ctrl_matrix
    file 'hf.bw' into bigWig_hf_matrix

    """
    export LD_LIBRARY_PATH=/usr/local/lib
    wiggletools mean ${bigwig_ctrl} | wigToBigWig stdin ${chrom_sizes} ctrl.bw
    wiggletools mean ${bigwig_hf} | wigToBigWig stdin ${chrom_sizes} hf.bw
    """
}

//return

// Parametrization for fast running, debugging
/*
before_start_length = 3000
body_length = 5000
after_end_length = 3000
*/
// Parametrization for production
before_start_length = 21600
body_length = 43200
after_end_length = 21600

// Order bigwig by group to show all mice of the same group together in the heatmap
bigWig_matrix.into { bigWig_matrix_c; bigWig_matrix_h }
bigWig_matrix_ctrl = bigWig_matrix_c
                        .filter( ~/.*[1,3,5,7,9].bw$/ )

bigWig_matrix_hf = bigWig_matrix_h
                        .filter( ~/.*[2,4,6,8,0].bw$/ )

/*
 * Creates deeptools matrix that will be used to generate a heatmap of the circadian rhythm by plotHeatmap
 */
process deep_tools_matrix {

    input:
    file dark_habituation from bed_dark_habituation
    file dark_development from bed_dark_development
    file bigwig_file_ctrl from bigWig_matrix_ctrl.toSortedList{ it.name.replace(".bw", "").toInteger() }
    file bigwig_file_hf from bigWig_matrix_hf.toSortedList{ it.name.replace(".bw", "").toInteger() }

    output:
    file 'matrix.mat.gz' into matrix_heatmap, matrix_profile

    """
    computeMatrix scale-regions -S ${bigwig_file_ctrl} ${bigwig_file_hf}\
                                -R ${dark_habituation} ${dark_development} \
                                --beforeRegionStartLength ${before_start_length} \
                                --regionBodyLength ${body_length} \
                                --afterRegionStartLength ${after_end_length} \
                                --skipZeros -out matrix.mat.gz
    """
}

/*
 *
 */
process deep_tools_matrix_groups {

    input:
    file dark_habituation from bed_dark_habituation_groups
    file dark_development from bed_dark_development_groups
    file bigwig_file_ctrl from bigWig_ctrl_matrix
    file bigwig_file_hf from bigWig_hf_matrix

    output:
    file 'matrix.mat.gz' into matrix_heatmap_group, matrix_profile_group

    """
    computeMatrix scale-regions -S ${bigwig_file_ctrl} ${bigwig_file_hf}\
                                -R ${dark_habituation} ${dark_development} \
                                --beforeRegionStartLength ${before_start_length} \
                                --regionBodyLength  ${body_length} \
                                --afterRegionStartLength ${after_end_length} \
                                --skipZeros -out matrix.mat.gz
    """
}

/*
 * Plots a heatmap comparing the circadian profile of the mice by group using deeptools
 */
process deep_tools_heatmap {

    publishDir "${params.output_res}/deeptools/", mode: 'copy', overwrite: 'true'

    input:
    file matrix from matrix_heatmap

    output:
    file "*.${image_format_deeptools}" into heatmap_fig

    """
    plotHeatmap -m ${matrix} \
                -out heatmap_actogram_like".${image_format_deeptools}" \
                -z Habituation Development \
                -T "Feeding behavior over 24 hours" \
                --startLabel APS \
                --endLabel RPS \
                --xAxisLabel "Time (s)" \
                --yAxisLabel "Food intake (g)" \
                --sortRegions no \
                --plotFileFormat ${image_format_deeptools} \
                --colorMap YlGnBu

    """
}

/*
 *
 */
process deep_tools_heatmap_by_groups {

    publishDir "${params.output_res}/deeptools/", mode: 'copy', overwrite: 'true'

    input:
    file matrix from matrix_heatmap_group

    output:
    file "*.${image_format_deeptools}" into heatmap_fig_gr

    """
    plotHeatmap -m ${matrix} \
                -out heatmap_actogram_like_groups".${image_format_deeptools}" \
                -z Habituation Development \
                -T "Feeding behavior over 24 hours" \
                --samplesLabel "Control mice" "High-Fat mice" \
                --startLabel APS \
                --endLabel RPS \
                --xAxisLabel "Time (s)" \
                --yAxisLabel "Food intake (g)" \
                --sortRegions no \
                --plotFileFormat ${image_format_deeptools} \
                --heatmapWidth 8 \
                --heatmapHeight 18 \
                --colorMap YlGnBu
    """
}

/*
 * Creates the circadian profile of the mice by group using deeptools
 */
process deep_tools_profile {

    publishDir "${params.output_res}/deeptools/", mode: 'copy', overwrite: 'true'

    input:
    file matrix from matrix_profile

    output:
    file "*.${image_format_deeptools}" into profile_fig

    """
    plotProfile -m ${matrix} \
                -out profile".${image_format_deeptools}" \
                --startLabel APS \
                --endLabel RPS \
                --yAxisLabel "Food intake (g)" \
                --plotTitle "Feeding behavior over 24 hours" \
                --regionsLabel Habituation Development \
                --plotFileFormat ${image_format_deeptools}
                #--xAxisLabel "Time (s)" \
    """
}

/*
 * Binning of meals based on meal duration, creates table files for chromHMM binarization and plots the distribution
 */
process bin {

    publishDir "${params.output_res}/chromHMM", mode: 'copy', pattern: "*.${image_format}", overwrite: 'true'

    input:
    file (dir_bed_feeding) from dir_bed_to_bin
    file 'cellmarkfiletable' from cell_mark_file_tbl

    output:
    file 'bed_binned' into dir_bed_binned
    file '*.binned' into cell_mark_file_tbl_binned
    file "meal_length_distro_binned.${image_format}" into plot_distro_binned

    """
    distro_meals_to_bin.R --path_bed_files=${dir_bed_feeding} \
                          --n_bins=${n_bins} \
                          --image_format=${image_format} > bins.txt

    for file_bed in ${dir_bed_feeding}/*.bed
    do
        bin_length_by_sliding_win.py -b \${file_bed} -bins "\$(< bins.txt)"
    done

    mkdir bed_binned
    mv *.bed bed_binned

    bins_string="\$(tr -d "\n\r" < bins.txt)"
    IFS=' ' read -r -a bins_ary <<< \$bins_string

    length_ary=\${#bins_ary[@]}
    i_last=\$((length_ary-1))
    i_for=\$((length_ary-2))

    awk -v bin_0=\${bins_ary[0]} \
        '{print \$1"\t0_"bin_0"_"\$2"\t0_"bin_0"_"\$3}' cellmarkfiletable > "${cellmarkfiletable}.binned"

    for index in \$(seq 0 \$i_for); do
        next_i=\$((index+1))
        awk -v bin_1=\${bins_ary[index]} -v bin_2=\${bins_ary[next_i]} \
            '{print \$1"\t"bin_1"_"bin_2"_"\$2"\t"bin_1"_"bin_2"_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    done

    awk -v  bin_l=\${bins_ary[i_last]} \
        '{print \$1"\t"bin_l"_"\$2"\t"bin_l"_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    """
}

/*
 * Binarizes feeding bed files (inputs 0 or 1 for each of the bins) using chromHMM
 */
process binarize {

    input:
    file chrom_sizes from chrom_sizes_chromHMM_binarize
    file dir_bed_binned from dir_bed_binned
    file 'cellmarkfiletable_binned' from cell_mark_file_tbl_binned

    output:
    file 'output_dir' into output_dir_binarized

    """
    mkdir output_dir

    java -mx4000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 \
                                                          -peaks ${chrom_sizes} \
                                                          ${dir_bed_binned} \
                                                          ${cellmarkfiletable_binned} \
                                                          output_dir
    """
}

/*
 * Learns the HMM model using chromHMM using the whole behavioral trajectory
 */
process HMM_model_learn {

    publishDir "${params.output_res}/chromHMM", mode: 'copy', pattern: "!(*.bed)    ", overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_l
    file 'input_binarized' from output_dir_binarized

    output:
    file 'output_learn/*dense*.bed' into HMM_model_ANNOTATED_STATES
    //file 'output_learn/*.*' into HMM_full_results
    file 'output_learn/emissions_*.*' into emission_results
    file 'output_learn/transitions_*.*' into transitions_results

    file '*.bed' into segmentation_bed_to_plot
    file '*dense.bed' into segmentation_bed_to_fraction
    file 'colormappingfile' into tbl_states_color

    """
    mkdir output_learn

    ## blue 1 active
    ## echo -e "1\t0,0,255" > colormappingfile
    ## red 2 resting
    ## echo -e "2\t255,0,0" >> colormappingfile
    ## yellow 3 snacking
    ## echo -e "3\t255,255,0" >> colormappingfile

    ## Using color blind palette
    # blue 1 active
    echo -e "1\t0,114,178" > colormappingfile
    # yellow 2 snacking
    echo -e "2\t240,228,66" >> colormappingfile
    # vermilion 3 resting
    echo -e "3\t213,94,0" >> colormappingfile
    # bluish green 4
    echo -e "4\t0,158,115" >> colormappingfile

    # java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes}  -printstatebyline test_feeding/output/outputdir test_feeding/output/outputdir_learn ${n_states} test_feeding/input/chrom.sizes
    java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 \
                                                         -l ${chrom_sizes} \
                                                         input_binarized output_learn \
                                                         ${n_states} ${chrom_sizes}

    for dense_file in output_learn/*segments*.bed
    do
        filename=\$(basename -- "\$dense_file")
        filename="\${filename%.*}"
        mice_id=\$(echo \$filename | cut -f2 -d_)

        java -mx4000M -jar /ChromHMM/ChromHMM.jar MakeBrowserFiles -c colormappingfile \
                                                                   \${dense_file} \
                                                                   \${mice_id} \
                                                                   \${filename}
    done
    """
}

/*
 * Uses gviz to plot the segmentation obtained with chromHMM
 */
process plot_HMM_states {

    publishDir "${params.output_res}/chromHMM/", mode: 'copy', overwrite: 'true'

    input:
    file 'output_learn/*' from segmentation_bed_to_plot.collect()
    file exp_phases_bed from exp_phases_bed_to_hmm
    file phases_bed from days_bed_gviz_hmm

    output:
    file "segmentation_HMM.${image_format}" into plot_HMM_segmentation

    """
    plot_HMM_segmentation.R --path_bed_files=output_learn \
                            --ini_time=0 \
                            --path_to_phases_file=${phases_bed} \
                            --path_to_exp_phases_file=${exp_phases_bed} \
                            --image_format=${image_format}
                            # --end_time=1814400 \
    """
}

/*
 * Calculates which is the fraction of time expend in each of the states during the behavioral trajectory
 */
process states_fraction {

    input:
    file (file_bed) from segmentation_bed_to_fraction.flatten()

    file exp_phases from exp_phases_bed_to_fraction.first()
    file dark_whole_experiment from whole_experiment_dark.first()
    file light_whole_experiment from whole_experiment_light.first()

    output:
    file '*.cov' into fraction_state_by_phase

    """
    fraction_states.py -b ${file_bed} -p ${dark_whole_experiment} -n ${n_states} -t "dark" # -c ${chrom_sizes}
    fraction_states.py -b ${file_bed} -p ${light_whole_experiment} -n ${n_states} -t "light" # -c ${chrom_sizes}
    """
}


/*
 * Compares the time expend in each state during habituation vs development in a plot
 */
process states_fraction_plots {

    publishDir "${params.output_res}/chromHMM/", mode: 'copy', overwrite: 'true'

    input:
    file './*' from fraction_state_by_phase.flatten().collect()
    file 'tbl_states_color' from tbl_states_color

    output:
    file "states_fraction.${image_format}" into plot_states_fraction

    """
    cat *.cov > all_bed_cov.txt
    states_fraction_time_plots.R --path_to_tbl=all_bed_cov.txt \
                                 --path_to_tbl_col=${tbl_states_color} \
                                 --image_format=${image_format}
    """
}
