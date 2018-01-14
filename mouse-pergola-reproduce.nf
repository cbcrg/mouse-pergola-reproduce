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

log.info "Mouse - Pergola - Reproduce  -  version 0.1.1"
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
  --exp_info='small_data/mappings/exp_info_small.txt' \
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
process stats_by_phase {

  	input:
  	file file_preferences from mice_files_preference
  	file mapping_file
    file mapping_file_phase
    file exp_phases

  	output:
  	//file 'stats_by_phase' into results_stats_by_phase
  	file 'behaviors_by_week' into d_behaviors_by_week
  	stdout into max_time
    file 'exp_phases' into exp_phases_bed_to_wr, exp_phases_bed_to_wr2
    file 'exp_phases_sushi' into exp_phases_bed_sushi, exp_phases_bed_gviz
	file 'stats_by_phase/phases_dark.bed' into exp_circadian_phases_sushi, exp_circadian_phases_gviz
    // stdout into test

  	"""
  	mice_stats_by_week.py -f "${file_preferences}"/intake*.csv -m ${mapping_file} -s "sum" -b feeding -p ${exp_phases} -mp ${mapping_file_phase}
  	mkdir stats_by_phase
  	mkdir behaviors_by_week
  	cp exp_phases.bed exp_phases
  	tail +2 exp_phases > exp_phases_sushi
  	mv *.bed  stats_by_phase/
  	mv *.tbl  behaviors_by_week/
  	"""
}


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

    if ( map_id_group.get("ctrl").contains(id.toInteger()) && food == "sc" )
    return "igv/1_ctrl_food_sc/${file}"

    if ( map_id_group.get("hf").contains(id.toInteger()) && food == "fat" )
    return  "igv/2_hf_food_fat/${file}"

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
  	for f in {1,3,5,7,11,13,15,17}
  	do
  	    mkdir -p work_dir

  	    files=( tr_"\$f"_* )

  	    if (( \${#files[@]} )); then
  	        cd work_dir
  	        track_int=`ls ../"tr_"\$f"_"*`
  	        mv \${track_int} track_int
  	        echo -e "food_sc\tblack" > dict_color
  	        echo -e "food_fat\tyellow" >> dict_color
  	        pergola_rules.py -i track_int -m ../${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color' -dl food_sc food_fat -d all
  	        in_f_sc=`ls tr_chr1*food_sc.bed`
            mv "\$in_f_sc" "`echo \$in_f_sc | sed s/chr1/\${f}/`"
  	        cd ..
  	        mv work_dir/tr*.bed ./
  	        mv work_dir/*.fa ./
        fi
  	done

    ## hf
  	for f in {2,4,8,10,12,14,16,18}
  	do
  	    mkdir -p work_dir

  	    files=( tr_"\$f"_* )

  	    if (( \${#files[@]} )); then
            cd work_dir
            track_int=`ls ../"tr_"\$f"_"*`
  	        mv \${track_int} track_int
            echo -e "food_sc\torange" > dict_color
  	        echo -e "food_fat\tblue" >> dict_color
            pergola_rules.py -i track_int -m ../${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color' -dl food_sc food_fat -d all
            # in_f_fat=`ls tr_chr1*food_fat.bed`
            in_f_fat=`ls tr_chr1*food_fat*.bed`
            mv "\$in_f_fat" "`echo \$in_f_fat | sed s/chr1/\${f}/`"
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

longest_phases_dark = phases_night
                        .max { it.size() }

longest_phases_dark.subscribe {
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
  	file 'tr*food*.bedGraph' into bedGraph_out, bedGraph_out_shiny_p, bedGraph_out_gviz, bedGraph_out_sushi
  	//file 'tr*{water,sac}*.bedGraph' into bedGraph_out_drink

  	"""
  	pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt -e -dl food_sc food_fat -d all
  	"""
}

bedGraph_out_shiny_p.flatten().subscribe {
    it.copyTo( result_dir_shiny_p.resolve ( ) )
}

/*
 * Plots preference of mice for the high-fat food in the different experimental phases
 */
/*
process plot_preference {

    publishDir "${params.output_res}/preference/", mode: 'copy', overwrite: 'true'

    input:
    file stats_by_phase from results_stats_by_phase

    output:
    file "*.${image_format}" into plot_preference

  	"""
    plot_preference.R --stat="sum" \
        --path2files=${stats_by_phase} \
        --image_format=${image_format}
  	"""
}
*/

exp_phases_bed_to_wr.subscribe {
    it.copyTo( result_dir_IGV.resolve ( 'exp_phases.bed' ) )
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