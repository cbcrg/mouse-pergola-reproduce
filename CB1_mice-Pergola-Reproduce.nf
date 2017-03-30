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
params.output      = "results/"

log.info "CB1_mice - Pergola - Reproduce  -  version 0.1"
log.info "====================================="
log.info "mice recordings        : ${params.recordings}"
log.info "mappings               : ${params.mappings}"
log.info "output                 : ${params.output}"
log.info "\n"

// Example command to run the script
/*
nextflow run CB1_mice-Pergola-Reproduce.nf \
  --recordings='small_data/mice_recordings/*.csv' \
  --mappings='small_data/mappings/b2p.txt' \
  -with-docker
*/
    
/*
 * Input parameters validation
 */
mapping_file = file(params.mappings)
mapping_file_bG = file(params.mappings)

/*
 * Input files validation
 */
if( !mapping_file.exists() ) exit 1, "Missing mapping file: ${mapping_file}"

/*
 * Create a channel for mice recordings 
 */
Channel
    .fromPath( params.recordings )
    .ifEmpty { error "Cannot find any CSV file with mice data" }
    .set { mice_files }

mice_files.into { mice_files_bed; mice_files_bedGraph }

process convert_bed {

  	input:
  	set file ('batch') from mice_files_bed
  	file mapping_file
  	
  	output:   	
  	set 'tr*food*.bed' into bed_out
  	set 'phases_dark.bed' into phases_dark
  	
  	"""  
  	pergola_rules.py -i ${batch} -m  ${mapping_file} -f bed -nt -e  	
  	"""
}

process convert_bedGraph {
  
  	input:
  	set file ('batch_bg') from mice_files_bedGraph
  	file mapping_file_bG
  	
  	output:   	
  	set 'tr*food*.bedGraph' into bedGraph_out
  	
  	"""  
  	# pergola_rules.py -i ${batch_bg} -m  ${mapping_file_bG} -f bedGraph -w 3600 -nt -e  
  	pergola_rules.py -i ${batch_bg} -m  ${mapping_file_bG} -f bedGraph -w 1800 -nt -e 	
  	"""
}

