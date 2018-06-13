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
params.phases_long    = "$baseDir/small_data/phases/exp_phases_whole_exp.csv"
params.mappings_phase = "$baseDir/small_data/mappings/f2g.txt"
params.exp_info       = "$baseDir/small_data/mappings/exp_info.txt"
params.output         = "files/"
params.image_format   = "tiff"

log.info "Mouse - Pergola - Reproduce  -  version 0.2.0"
log.info "====================================="
log.info "mice recordings          : ${params.recordings}"
log.info "mappings                 : ${params.mappings}"
log.info "mappings bed             : ${params.mappings_bed}"
log.info "experimental phases      : ${params.phases}"
log.info "experimental phases long : ${params.phases_long}"
log.info "mappings phases          : ${params.mappings_phase}"
log.info "experimental info        : ${params.exp_info}"
log.info "output                   : ${params.output}"
log.info "image format             : ${params.image_format}"
log.info "chromHMM ct config table : ${params.tbl_chromHMM_ctrl}"
log.info "chromHMM hf config table : ${params.tbl_chromHMM_hf}"
log.info "\n"

// Example command to run the script
/*
nextflow run mouse-pergola-reproduce.nf \
  --recordings='small_data/mouse_recordings/' \
  --mappings='small_data/mappings/b2p.txt' \
  --mappings_bed='small_data/mappings/bed2pergola.txt' \
  --phases='small_data/phases/exp_phases.csv' \
  --phases_long='small_data/phases/exp_phases_whole_exp.csv' \
  --mappings_phase='small_data/mappings/f2g.txt' \
  --exp_info='small_data/mappings/exp_info.txt' \
  --image_format='png' \
  --tbl_chromHMM_ctrl="small_data/chromHMM_files/cellmarkfiletable_ctrl" \
  --tbl_chromHMM_hf="small_data/chromHMM_files/cellmarkfiletable_hf" \
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
exp_phases_long = file(params.phases_long)
exp_info = file(params.exp_info)

cell_mark_file_tbl_ctrl = Channel.fromPath(params.tbl_chromHMM_ctrl)
cell_mark_file_tbl_hf = Channel.fromPath(params.tbl_chromHMM_hf)

/*
 * Input files validation
 */
if( !mapping_file.exists() ) exit 1, "Missing mapping file: ${mapping_file}"
if( !mapping_file_phase.exists() ) exit 1, "Missing mapping phases file: ${mapping_file_phase}"
if( !exp_phases.exists() ) exit 1, "Missing phases file: ${exp_phases}"
if( !exp_phases_long.exists() ) exit 1, "Missing long phases file: ${exp_phases_long}"
if( !exp_info.exists() ) exit 1, "Missing experimental info file: ${exp_info}"
//if( !cell_mark_file_tbl_ctrl.first().exists() ) exit 1, "Missing configuration file for chromHMM ctrl: ${cell_mark_file_tbl_ctrl}"
//if( !cell_mark_file_tbl_hf.first().exists() ) exit 1, "Missing configuration file for chromHMM hf: ${cell_mark_file_tbl_hf}"

cell_mark_file_tbl_ctrl_labelled = cell_mark_file_tbl_ctrl.map {
                                        [ it, "ctrl" ]
                                   }
cell_mark_file_tbl_hf_labelled = cell_mark_file_tbl_hf.map {
                                        [ it, "hf" ]
                                   }

cell_mark_file_tbl = cell_mark_file_tbl_ctrl_labelled.mix ( cell_mark_file_tbl_hf_labelled )

/*
 * Read image format
 */
image_format = "${params.image_format}"

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
    file exp_phases_long

  	output:
  	file 'behaviors_by_week' into d_behaviors_by_week
  	stdout into max_time
    file 'exp_phases' into exp_phases_bed_to_wr, exp_phases_bed_to_wr2
    file 'exp_phases_sushi' into exp_phases_bed_sushi, exp_phases_bed_gviz

	file 'stats_by_phase/phases_dark.bed' into exp_circadian_phases_sushi, exp_circadian_phases_gviz, days_bed_igv, days_bed_shiny, days_bed_deepTools
    file 'Habituation_dark_long.bed' into bed_dark_habituation
    file 'Development_dark_long.bed' into bed_dark_development

  	"""
  	mice_stats_by_week.py -f "${file_preferences}"/intake*.csv -m ${mapping_file} -s "sum" -b feeding -p ${exp_phases} \
  	                      -pl ${exp_phases_long} -mp ${mapping_file_phase}

  	mkdir stats_by_phase
  	mkdir behaviors_by_week
  	cp exp_phases.bed exp_phases
  	tail +2 exp_phases > exp_phases_sushi
  	mv *.bed  stats_by_phase/
  	mv stats_by_phase/Development_dark_long.bed ./
  	mv stats_by_phase/Habituation_dark_long.bed ./
  	mv *.tbl  behaviors_by_week/
  	"""
}

/*
 * Creates a heatmap that compares mice feeding behavior of high-fat mice with their controls
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
  	file 'dir_bed_ctrl' into dir_bed_to_chromHMM_ctrl
  	file 'dir_bed_hf' into dir_bed_to_chromHMM_hf
  	// file 'tr*{water,sac}*.bed' into bed_out_drink
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


  	"""
}

dir_bed_to_chromHMM_ctrl_labelled = dir_bed_to_chromHMM_ctrl.map {
                                        [ it, "ctrl" ]
                                    }

dir_bed_to_chromHMM_hf_labelled = dir_bed_to_chromHMM_hf.map {
                                    [ it, "hf" ]
                                  }

dir_bed_to_chromHMM = dir_bed_to_chromHMM_ctrl_labelled.mix ( dir_bed_to_chromHMM_hf_labelled )

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
  	file 'tr*food*.bedGraph' into bedGraph_out, bedGraph_out_shiny_p, bedGraph_out_gviz, bedGraph_out_sushi, bedGraph_out_bigwig
  	file 'chrom.sizes' into chrom_sizes, chrom_sizes_chromHMM_bin, chrom_sizes_chromHMM_binarize, chrom_sizes_chromHMM_l
  	//file 'tr*{water,sac}*.bedGraph' into bedGraph_out_drink
    //stdout into len_experiment_days
    stdout into len_experiment_weeks

  	"""
  	pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt -e -dl food_sc food_fat -d all
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

//bedGraph_out_bigwig.flatten().println()

bedGraph_out_bigwig_short_name = bedGraph_out_bigwig.flatten().map {
    def content = it
    def name = it.baseName.replaceAll('tr_','').replaceAll('_dt_food_fat_food_sc','').replaceAll('_dt_food_sc','')
    [ content, name ]
}


/*
 * For each pair get the correlation using bigwig
 */
process bedgraph_to_bigWig {
    //container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

  	input:
  	set file (bedgr_file), val (name) from bedGraph_out_bigwig_short_name
    file chrom_sizes from chrom_sizes.first()

    output:
    file '*.bw' into bigWig_matrix, bigWig_multiBigwigSummary

    """
  	head -n -2 ${bedgr_file} > ${name}.trimmed
    bedGraphToBigWig ${name}.trimmed ${chrom_sizes} ${name}".bw"
  	"""
}

process deep_tools_matrix {
    //container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:
    //file day_phases_dir from days_bed_deepTools
    file dark_habituation from bed_dark_habituation
    file dark_development from bed_dark_development
    file bigwig_file from bigWig_matrix.toSortedList{ it.name.replace(".bw", "").toInteger() }

    output:
    file 'matrix.mat.gz' into matrix_heatmap, matrix_profile

    """
    computeMatrix scale-regions -S ${bigwig_file} \
                                -R  ${dark_habituation} ${dark_development} \
                                --beforeRegionStartLength 3000 \
                                --regionBodyLength  5000 \
                                --afterRegionStartLength 3000 \
                                --skipZeros -out matrix.mat.gz



    """
}

/*
nice plot
computeMatrix scale-regions -S ${bigwig_file} \
                                -R phases_dark.bed \
                                -R  ${dark_habituation} ${dark_development} \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -out matrix.mat.gz
whole plot
                              --beforeRegionStartLength 21600 \
                              --regionBodyLength  43200 \
                              --afterRegionStartLength 21600 \
*/

//file '*.png' into heatmap_deepTools

process deep_tools_heatmap {
    //container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:
    file matrix from matrix_heatmap

    output:
    file "*.${image_format}" into heatmap_fig

    """
    plotHeatmap -m ${matrix} \
                -out heatmap_actogram_like".${image_format}" \
                -z Habituation Development \
                -T "Feeding behavior over 24 hours" \
                --startLabel APS \
                --endLabel RPS \
                --xAxisLabel "Time (s)" \
                --plotFileFormat ${image_format} #\
                #--sortRegions no
                # --kmeans 2


    """
}

process deep_tools_profile {
    //container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:
    file matrix from matrix_profile

    output:
    file "*.${image_format}" into profile_fig

    """
    plotProfile -m ${matrix} \
                -out profile".${image_format}" \
                --plotFileFormat ${image_format} \
                --plotTitle "Test data profile"


    """
}

/*
process deep_tools_bigWigSummary {
    //container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:
    file bigwig_file from bigWig_multiBigwigSummary.toSortedList{ it.name.replace(".bw", "").toInteger() }

    output:
    file 'results.npz' into multiBigwigSummary
    file "PCA.${image_format}" into pca_fig

    """
    multiBigwigSummary bins -b ${bigwig_file} \
                            -o results.npz #\
                            #-bs 1800 #\
                            #--smartLabels # label is file without extension

    plotPCA -in results.npz \
            -o PCA".${image_format}" \
            -T "PCA" \
            --plotFileFormat ${image_format} #\
            #--colors #999999 #e69f00 #999999 #e69f00 #999999 #999999 #e69f00 #999999 #e69f00 #999999 #e69f00 #999999 #e69f00 #999999 #e69f00 #999999 #e69f00
    """
}
*/

/*
process deep_tools_pca {
    //container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:
    file bigwig_file from bigWig_multiBigwigSummary.toSortedList{ it.name.replace(".bw", "").toInteger() }

    output:
    file 'results.npz' into multiBigwigSummary

    """
    multiBigwigSummary -b $bigwig_file -o results.npz


    """
}
*/

//step = 3600 * 24
step = 604800
window = 604800
//len_days = len_experiment_days.getVal().toInteger()
len_weeks = len_experiment_weeks.getVal().toInteger()
start_end_win = Channel
                    .from( 0..len_weeks )
                    //.from( 0..len_days )
                    .map { [ it, 1 + (it * step),  + window + (it * step) ] }

//start_end_win.println()

// para ordenar actograms
/*
Channel
    .from( 1, 2, 3, 4, 5 )
    .filter { it % 2 == 1 }
*/


//bins="30 120"
bins="30"
bins_tbl = bins

dir_bed_to_chromHMM_win = dir_bed_to_chromHMM.spread (start_end_win)
//dir_bed_to_chromHMM_win.into { culo; dir_bed_to_chromHMM_win } //del
//culo.println() //del

/*
 * Intersects the data with the given start and end and produces the files containing the meals corresponding to a
 * given bin
 */
process bin {

    publishDir "results/", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_bin.first()
    set file (dir_bed_feeding), val (group), val (index), val (start), val (end) from dir_bed_to_chromHMM_win
    val bins

    output:
    set index, group, 'output_bed_binned' into output_dir_binned

    """
    rm -f [0-9]*.bed

    for file_bed in ${dir_bed_feeding}/*.bed
    do
        # bin_length_by_sliding_win.py -b \${file_bed} -ct 1 604800 -bins 30 120
        bin_length_by_sliding_win.py -b \${file_bed} -ct ${start} ${end} -bins ${bins}
    done
    echo ${index}
    mkdir output_bed_binned
    mv *.bed output_bed_binned
    """
}

process create_chromHMM_table {

    publishDir "results/", mode: 'copy', overwrite: 'true'

    input:
    set file ('cellmarkfiletable'), val (group) from cell_mark_file_tbl
    val bins_tbl

    output:
    set '*.binned', group into cell_mark_file_tbl_binned

    """
    bins_ary=(${bins})
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

output_dir_binned_and_cell_mark_tbl = output_dir_binned.spread ( cell_mark_file_tbl_binned )
                                        .filter { it[1] == it[4] }
                                        .map { [ it[0], it[1], it[2], it[3] ] }


/*
 * chromHMM binarizes feeding bed files
 */
process binarize {
    publishDir "results/", mode: 'copy', overwrite: 'true'
    memory 1.GB

    input:
    file chrom_sizes from chrom_sizes_chromHMM_binarize.first()
    //set val (index), val (group), file (dir_bed_feeding) from output_dir_binned
    //set file ('cellmarkfiletable'), val (group) from cell_mark_file_tbl_binned
    set val (index), val (group), file (dir_bed_feeding), file ('cellmarkfiletable') from output_dir_binned_and_cell_mark_tbl

    output:
    set index, group, 'output_dir' into output_dir_binarized

    """
    cat ${cellmarkfiletable} > ${cellmarkfiletable}".filtered"

    mkdir output_dir
    # java -mx1000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 -peaks ${chrom_sizes} ${dir_bed_feeding} ${cellmarkfiletable} output_dir
    java -mx1000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 -peaks ${chrom_sizes} ${dir_bed_feeding} ${cellmarkfiletable}".filtered" output_dir
    """
}

/*
if [ "$index" -eq "2" ]
    then
        cat ${cellmarkfiletable} | grep -v "tr_18" | grep -v "tr_10" | grep -v "tr_12" > ${cellmarkfiletable}".filtered"
    fi

    if [ "$index" -eq "4" ]
    then
        # cat ${cellmarkfiletable} | grep -v "tr_7" | grep -v "tr_17" > ${cellmarkfiletable}".filtered" # ok
        cat ${cellmarkfiletable} | grep -v "tr_7" > ${cellmarkfiletable}".filtered" # con esto medio bien pero hay pico
    fi

    if [ "$index" -eq "6" ]
    then
        cat ${cellmarkfiletable} |  grep -v "tr_3" |  grep -v "tr_17" > ${cellmarkfiletable}".filtered"
        #  grep -v "tr_1_"
    fi

*/

/*
n_states = 2

process HMM_model_learn {

    publishDir "${params.output_res}/chromHMM", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_l.first()
    set val (index), val (group), file ('input_binarized') from output_dir_binarized

    output:
    set group, "output_learn_${group}/*dense*.bed" into HMM_model_ANNOTATED_STATES
    set group, index, "output_learn_${group}/*.*" into HMM_full_results
    file 'model.tbl' into model

    """
    mkdir output_learn_${group}
    # java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes}  -printstatebyline test_feeding/output/outputdir test_feeding/output/outputdir_learn ${n_states} test_feeding/input/chrom.sizes
    java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes} input_binarized output_learn_${group} ${n_states} ${chrom_sizes}

    # cat "emissions_"${n_states}".txt" | head -1 | awk -v index_v=${index} '{index_v"\t"\$0}' > emissions.tbl

    tail +2 output_learn_${group}"/model_"${n_states}".txt" | awk -v index_v=${index} \
                                                                  -v group=${group} '{print index_v"\t"group"\t"\$0}' >> model.tbl

    ## change color of annotations
    for dense_file in output_learn_${group}/*segments*.bed
    do
        filename=\$(basename -- "\$dense_file")
        filename="\${filename%.*}"
        mice_id=\$(echo \$filename | cut -f2 -d_)

        java -mx4000M -jar /ChromHMM/ChromHMM.jar MakeBrowserFiles -c colormappingfile \${dense_file} \${mice_id} \${filename} ${chrom_sizes}
    done
    """
}
*/

ctrl = Channel
    .from( 1,3,5,7,9,11,13,15,17 )

hf = Channel
    .from( 2,4,8,10,12,14,16,18 )

output_dir_binarized.into { output_dir_binarized; output_dir_binarized_to_ctrl;  output_dir_binarized_to_hf }

output_dir_binarized_ctrl = output_dir_binarized_to_ctrl.filter { it[1] == 'ctrl' }
output_dir_binarized_hf = output_dir_binarized_to_hf.filter { it[1] == 'hf' }

output_dir_to_binarized_ctrl = output_dir_binarized_ctrl.spread (ctrl)
//output_dir_to_binarized_ctrl.println()

output_dir_to_binarized_hf = output_dir_binarized_hf.spread (hf)
//output_dir_to_binarized_hf.println()

//dir_bed_to_chromHMM = dir_bed_to_chromHMM_ctrl_labelled.mix ( dir_bed_to_chromHMM_hf_labelled )


//este
output_dir_to_binarized_all = output_dir_to_binarized_ctrl.mix ( output_dir_to_binarized_hf )

//dir_bed_to_chromHMM_win.into { culo; dir_bed_to_chromHMM_win } //del
//output_dir_to_binarized_ctrl.println() //del
//output_dir_to_binarized_hf.println() //del

//output_dir_to_binarized_all.println() //del

n_states = 2

process HMM_model_learn {
    tag { mice_id + "_" + index }

    publishDir "${params.output_res}/chromHMM", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_l.first()
    set val (index), val (group), file ('input_binarized'), val (mice_id) from output_dir_to_binarized_all

    output:
    set group, "output_learn/*dense*.bed" into HMM_model_ANNOTATED_STATES
    set group, index, "output_learn/*.*" into HMM_full_results
    file 'model.tbl' into model

    """
    mkdir output_learn
    mkdir input_binarized_one_mouse

    if [ ! -f input_binarized/tr_${mice_id}_chr1_binary.txt ]; then
        touch output_learn/tr_${mice_id}dense_empty.bed
        touch model.tbl
        #echo -e "1\t$group\t${mice_id}\ttransitionprobs\t1\t1\tNA" > model1.tbl
        #echo -e "1\t$group\t${mice_id}\ttransitionprobs\t1\t2\tNA" >> model1.tbl
        #echo -e "1\t$group\t${mice_id}\ttransitionprobs\t2\t2\tNA" >> model1.tbl
        #echo -e "1\t$group\t${mice_id}\ttransitionprobs\t2\t1\tNA" >> model1.tbl

    else
        cp input_binarized/tr_${mice_id}_* input_binarized_one_mouse

        java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes} input_binarized_one_mouse output_learn ${n_states} ${chrom_sizes}

        tail +2 output_learn"/model_"${n_states}".txt" | awk -v index_v=${index} \
                                                             -v group=${group} \
                                                             -v mice_id=${mice_id} \
                                                             '{print index_v"\t"group"\t"mice_id"\t"\$0}' >> model.tbl
    fi

    ## change color of annotations
    # for dense_file in output_learn/*segments*.bed
    # do
    #   filename=\$(basename -- "\$dense_file")
    #   filename="\${filename%.*}"
    #   mice_id=\$(echo \$filename | cut -f2 -d_)
    #
    #   java -mx4000M -jar /ChromHMM/ChromHMM.jar MakeBrowserFiles -c colormappingfile \${dense_file} \${mice_id} \${filename} ${chrom_sizes}
    # done
    """
}


process subset_tbl_by_trans_emission {
    publishDir "${params.output_res}/chromHMM", mode: 'copy', overwrite: 'true'

    input:
    file 'model.tbl' from model.collectFile()

    output:
    file 'emission_prob.tbl' into emissions
    file 'transition_prob.tbl' into transitions
    file 'prob_init.tbl' into probs_init

    """
    tail +2 model.tbl | grep "emissionprobs" >> emission_prob.tbl
    tail +2 model.tbl | grep "probinit" >> prob_init.tbl
    tail +2 model.tbl | grep "transitionprobs" >> transition_prob.tbl
    """

}

// emissions.println()
// transitions.println()
/*
process plot_transistions {
    publishDir "${params.output_res}/results_segmentation", mode: 'copy', overwrite: 'true'

    input:
    file emissions_tbl from emissions
    file transitions_tbl from transitions

    output:

    file "transition.${image_format}" into transition_plot
    file "emission.${image_format}" into emission_plot

    """
    plot_em_trans_by_win.R --path2ems=${emissions_tbl} \
                           --path2trans=${transitions_tbl} \
                           --image_format=${image_format}
    """

}
*/
/*
# awk '{print \$1"\t""30_"\$2"\t""30_"\$3}' cellmarkfiletable > "${cellmarkfiletable}.binned"
# awk '{print \$1"\t""30_120"\$2"\t""30_120"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    # awk '{print \$1"\t""120_"\$2"\t""120_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
*/
/*
# bin_length_by_sliding_win.py -b \${file_bed} -ct ${start} ${end} -bins 30 120

    # pipe bedtools intersect with the respective window

    #for day in {1..last_day}

    # for day in length
        # crear un pequenyo script the python que lo haga, tiene que tener como input:
            # start and end of bed representing a week
            # bed file para intersectar
            # los bins que voy a usar [0,30], [30,120], [120]
        # generate bed of day
            echo -e "1\t \t\n"
        # primero intersectar o primero hacer el awk, comparar running time
        # primero intersecto luego el awk actua sobre un archivo mas corto y sera mas rapido
        #for file_bed in *.bed
        #do
        #    bedtools intersect -a \${file_bed}  -b bed_day_\${n} > bed_intersect

        #done




    for file_bed in *.bed
    do
    # esta parte quiza sea mas facil hacerla aqui
	    cat \${file_bed} | awk -F'\t' '\$3-\$2<=30 {print \$0}' > "30_\${file_bed}"
	    cat \${file_bed} | awk -F'\t' '\$3-\$2>30 && \$3-\$2<=120{print \$0}' > "30_120\${file_bed}"
	    cat \${file_bed} | awk -F'\t' '\$3-\$2>120 {print \$0}' > "120_\${file_bed}"
    done

    cd ..

    awk '{print \$1"\t""30_"\$2"\t""30_"\$3}' cellmarkfiletable > "${cellmarkfiletable}.binned"
    awk '{print \$1"\t""30_120"\$2"\t""30_120"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    awk '{print \$1"\t""120_"\$2"\t""120_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"

    mkdir output_dir

    java -mx4000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 -peaks ${chrom_sizes} ${dir_bed_feeding} "${cellmarkfiletable}.binned" output_dir
    """
}


/*
 * chromHMM binarizes feeding bed files
 */
/*
process bin_binarize {
    publishDir "results/", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_b
    file dir_bed_feeding from dir_bed_to_chromHMM
    file 'cellmarkfiletable' from cell_mark_file_tbl

    output:
    file 'output_dir' into output_dir_binarized

    """
    cd ${dir_bed_feeding}

    rm -f [0-9]*.bed

    for file_bed in *.bed
    do
	    cat \${file_bed} | awk -F'\t' '\$3-\$2<=30 {print \$0}' > "30_\${file_bed}"
	    cat \${file_bed} | awk -F'\t' '\$3-\$2>30 && \$3-\$2<=120{print \$0}' > "30_120\${file_bed}"
	    cat \${file_bed} | awk -F'\t' '\$3-\$2>120 {print \$0}' > "120_\${file_bed}"
    done

    cd ..

    awk '{print \$1"\t""30_"\$2"\t""30_"\$3}' cellmarkfiletable > "${cellmarkfiletable}.binned"
    awk '{print \$1"\t""30_120"\$2"\t""30_120"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    awk '{print \$1"\t""120_"\$2"\t""120_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"

    mkdir output_dir

    java -mx4000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 -peaks ${chrom_sizes} ${dir_bed_feeding} "${cellmarkfiletable}.binned" output_dir
    """
}
*/
/*
 * chromHMM binarizes feeding bed files
 */
/*
process binarize {
    publishDir "results/", mode: 'copy', overwrite: 'true'
    memory 4.GB

    input:
    file chrom_sizes from chrom_sizes_chromHMM_b
    file dir_bed_feeding from dir_bed_to_chromHMM
    file cellmarkfiletable from cell_mark_file_tbl

    output:
    file 'output_dir' into output_dir_binarized

    """
    mkdir output_dir
    java -mx4000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 -peaks ${chrom_sizes} ${dir_bed_feeding} ${cellmarkfiletable} output_dir
    """
}

process HMM_model_learn {

    publishDir "${params.output_res}/chromHMM", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_l
    file 'input_binarized' from output_dir_binarized

    output:
    file '*dense*.bed' into HMM_model_ANNOTATED_STATES
    file 'output_learn/*.*' into HMM_full_results

    """
    mkdir output_learn
    # java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes}  -printstatebyline test_feeding/output/outputdir test_feeding/output/outputdir_learn ${n_states} test_feeding/input/chrom.sizes
    java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes} input_binarized output_learn ${n_states} ${chrom_sizes}

    ## change color of annotations
    for dense_file in output_learn/*segments*.bed
    do
        filename=\$(basename -- "\$dense_file")
        filename="\${filename%.*}"
        mice_id=\$(echo \$filename | cut -f2 -d_)

        java -mx4000M -jar /ChromHMM/ChromHMM.jar MakeBrowserFiles -c colormappingfile \${dense_file} \${mice_id} \${filename}
    done
    """
}


*/
