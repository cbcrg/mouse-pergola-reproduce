#!/usr/bin/env nextflow

/*
 *  Copyright (c) 2014-2017, Centre for Genomic Regulation (CRG).
 *  Copyright (c) 2014-2017, Jose Espinosa-Carrasco and the respective authors.
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

params.recordings  = "$baseDir/small_data/mice_recordings/*.csv"
params.mappings    = "$baseDir/small_data/mappings/b2p.txt"   
params.output      = "files/"

log.info "CB1_mice - Pergola - Reproduce  -  version 0.1"
log.info "====================================="
log.info "mice recordings        : ${params.recordings}"
log.info "mappings               : ${params.mappings}"
log.info "experimental phases    : ${params.phases}"
log.info "mappings phases        : ${params.mappings_phase}"
log.info "output                 : ${params.output}"
log.info "\n"

// Example command to run the script
/*
nextflow run CB1_mice-Pergola-Reproduce.nf \
  --recordings='small_data/mice_recordings/' \
  --mappings='small_data/mappings/b2p.txt' \
  --phases='small_data/mice_recordings/exp_phases.csv' \
  --mappings_phase='small_data/mappings/f2g.txt' \
  -with-docker
*/
    
/*
 * Input parameters validation
 */
mapping_file = file(params.mappings)
mapping_file_bG = file(params.mappings)
mapping_file_phase = file(params.mappings_phase)

exp_phases = file(params.phases)

/*
 * Input files validation
 */
if( !mapping_file.exists() ) exit 1, "Missing mapping file: ${mapping_file}"
if( !mapping_file_phase.exists() ) exit 1, "Missing mapping phases file: ${mapping_file_phase}"
if( !exp_phases.exists() ) exit 1, "Missing phases file: ${exp_phases}"

/*
 * Create a channel for mice recordings 
 */
Channel
    .fromPath( "${params.recordings}intake*.csv" )
    //"${params.var_dir_test}/*.txt"
    .ifEmpty { error "Cannot find any CSV file with mice data" }
    .set { mice_files }

mice_files.into { mice_files_bed; mice_files_bedGraph }

/*
 * Create a channel for mice recordings
 */

Channel
    .fromPath( params.recordings )
    .set { mice_files_preference }

//mice_files_preference.println()

process convert_bed {

    publishDir = [path: {params.output}, mode: 'copy', overwrite: 'true']

  	input:
  	file ('batch') from mice_files_bed
  	file mapping_file
  	
  	output:   	
  	file 'tr*food*.bed' into bed_out
  	file 'phases_dark.bed' into phases_dark
  	
  	"""  
  	pergola_rules.py -i ${batch} -m  ${mapping_file} -f bed -nt -e
  	"""
}

process convert_bedGraph {

    publishDir = [path: {params.output}, mode: 'copy', overwrite: 'true']

  	input:
  	file ('batch_bg') from mice_files_bedGraph
  	file mapping_file_bG
  	
  	output:   	
  	file 'tr*food*.bedGraph' into bedGraph_out
  	
  	"""  
  	# pergola_rules.py -i ${batch_bg} -m  ${mapping_file_bG} -f bedGraph -w 3600 -nt -e  
  	pergola_rules.py -i ${batch_bg} -m  ${mapping_file_bG} -f bedGraph -w 1800 -nt -e 	
  	"""
}

process stats_by_phase {

  	input:
  	set file_preferences from mice_files_preference
  	file mapping_file
    file mapping_file_phase
    file exp_phases

  	output:
  	file 'stats_by_phase' into results_stats_by_phase
    //set 'results_score', tag_group into results_bed_score, results_bed_score_2

  	"""
  	mice_stats_by_phase.py -f "${file_preferences}"/intake*.csv -m ${mapping_file} -s "sum" -b feeding -p ${exp_phases} -mp ${mapping_file_phase}
  	mkdir stats_by_phase
  	mv *.bed  stats_by_phase/
  	"""
}

process plot_preference {

    publishDir = [path: "plot", mode: 'copy', overwrite: 'true']

    input:
    file stats_by_phase from results_stats_by_phase

    output:
    file '*.png' into plot_preference

  	"""
    plot_preference.R --stat="sum" --path2files=${stats_by_phase}
    #--path2plot="/Users/jespinosa/git/pergola-paper-reproduce/cb1_mice/bin/results/"
  	"""
}
