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
params.output         = "files/"
params.image_format   = "tiff"

log.info "Mouse - Pergola - Reproduce  -  version 0.2.0"
log.info "====================================="
log.info "mice recordings        : ${params.recordings}"
log.info "mappings               : ${params.mappings}"
log.info "mappings bed           : ${params.mappings_bed}"
log.info "experimental phases    : ${params.phases}"
log.info "mappings phases        : ${params.mappings_phase}"
log.info "experimental info      : ${params.exp_info}"
log.info "output                 : ${params.output}"
log.info "image format           : ${params.image_format}"
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
  --image_format='tiff' \
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

  	output:
  	file 'behaviors_by_week' into d_behaviors_by_week
  	stdout into max_time
    file 'exp_phases' into exp_phases_bed_to_wr, exp_phases_bed_to_wr2
    file 'exp_phases_sushi' into exp_phases_bed_sushi, exp_phases_bed_gviz
	file 'stats_by_phase/phases_dark.bed' into exp_circadian_phases_sushi, exp_circadian_phases_gviz, days_bed_igv, days_bed_shiny, days_bed_deepTools
    file 'Development_dark.bed' into bed_dark_development
    file 'Habituation_dark.bed' into bed_dark_habituation

  	"""
  	mice_stats_by_week.py -f "${file_preferences}"/intake*.csv -m ${mapping_file} -s "sum" -b feeding -p ${exp_phases} -mp ${mapping_file_phase}
  	mkdir stats_by_phase
  	mkdir behaviors_by_week
  	cp exp_phases.bed exp_phases
  	tail +2 exp_phases > exp_phases_sushi
  	mv *.bed  stats_by_phase/
  	mv stats_by_phase/Development_dark.bed ./
  	mv stats_by_phase/Habituation_dark.bed ./
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

    def map_id_group = [ "ctrl" : [1,3,5,7,11,13,15,17],
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
  	// file 'tr*{water,sac}*.bed' into bed_out_drink
  	file 'phases_light.bed' into phases_night
  	file '*.fa' into out_fasta

  	"""
  	pergola_rules.py -i ${batch} -m ${mapping_file} -f bed -nt -e -bl -dl food_sc food_fat -d all

    shopt -s nullglob

	## ctrl
  	for f in {1,3,,5,7,11,13,15,17}
  	do
  	    mkdir -p work_dir

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
  	        mv work_dir/tr*.bed ./
  	        mv work_dir/*.fa ./
        fi
  	done

    for f in {2,4,8,10,12,14,16,18}
  	do
  	    mkdir -p work_dir

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
  	        mv work_dir/tr*.bed ./
  	        mv work_dir/*.fa ./
        fi
  	done
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
  	file 'tr*food*.bedGraph' into bedGraph_out, bedGraph_out_shiny_p, bedGraph_out_gviz, bedGraph_out_sushi, bedGraph_out_bigwig
  	file 'chrom.sizes' into chrom_sizes
  	//file 'tr*{water,sac}*.bedGraph' into bedGraph_out_drink

  	"""
  	pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt -e -dl food_sc food_fat -d all
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
    container = '089a918d085e'
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
//bigWig_out.toSortedList{ it.name.replace(".bw", "").toInteger() }.println()

process deep_tools_matrix {
    container = '089a918d085e'
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
                              --beforeRegionStartLength 21600 \
                              --regionBodyLength  43200 \
                              --afterRegionStartLength 21600 \
                              --skipZeros -out matrix.mat.gz



    """
}

/*
nice plot
computeMatrix scale-regions -S ${bigwig_file} \
                                -R phases_dark.bed \
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
    container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:
    file matrix from matrix_heatmap

    output:
    file '*.png' into heatmap_fig

    """
    plotHeatmap -m ${matrix} \
                -out heatmap_actogram_like.png \
                --kmeans 2


    """
}

process deep_tools_profile {
    container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:
    file matrix from matrix_profile

    output:
    file '*.png' into profile_fig

    """
    plotProfile -m ${matrix} \
                -out profile.png \
                --plotTitle "Test data profile"


    """
}

process deep_tools_bigWigSummary {
    container = '089a918d085e'
    publishDir "results_new/", mode: 'copy', overwrite: 'true'

    input:

    file bigwig_file from bigWig_multiBigwigSummary.toSortedList{ it.name.replace(".bw", "").toInteger() }

    output:


    file 'results.npz' into multiBigwigSummary
    file 'PCA.png' into pca_fig

    """
    multiBigwigSummary bins -b $bigwig_file -o results.npz

    plotPCA -in results.npz \
            -o PCA.png \
            -T "PCA"
    """
}

/*
process deep_tools_pca {
    container = '089a918d085e'
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