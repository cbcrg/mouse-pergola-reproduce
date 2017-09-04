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
log.info "experimental info      : ${params.exp_info}"
log.info "output                 : ${params.output}"
log.info "\n"

// Example command to run the script
/*
nextflow run CB1_mice-Pergola-Reproduce.nf \
  --recordings='small_data/mice_recordings/' \
  --mappings='small_data/mappings/b2p.txt' \
  --mappings_bed='small_data/mappings/bed2pergola.txt' \
  --phases='small_data/mice_recordings/exp_phases.csv' \
  --mappings_phase='small_data/mappings/f2g.txt' \
  --exp_info='small_data/mappings/exp_info_small.txt' \
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

process stats_by_phase {

  	input:
  	file file_preferences from mice_files_preference
  	file mapping_file
    file mapping_file_phase
    file exp_phases

  	output:
  	file 'stats_by_phase' into results_stats_by_phase
  	stdout into max_time
    file 'exp_phases' into exp_phases_bed_gviz, exp_phases_bed_to_wr
    file 'exp_phases_sushi' into exp_phases_bed_sushi

  	"""
  	mice_stats_by_phase.py -f "${file_preferences}"/intake*.csv -m ${mapping_file} -s "sum" -b feeding -p ${exp_phases} -mp ${mapping_file_phase}
  	mkdir stats_by_phase
  	cp exp_phases.bed exp_phases
  	tail +2 exp_phases > exp_phases_sushi
  	mv *.bed  stats_by_phase/
  	"""
}

process convert_bed {

    publishDir = [path: {params.output}, mode: 'copy', overwrite: 'true']

  	input:
  	file ('batch') from mice_files_bed
  	file mapping_file
  	file mapping_bed_file
  	
  	output:   	
  	file 'tr*food*.bed' into bed_out
  	file 'tr*{water,sac}*.bed' into bed_out_drink
  	file 'phases_dark.bed' into phases_dark

  	"""
  	pergola_rules.py -i ${batch} -m ${mapping_file} -f bed -nt -e

    shopt -s nullglob

  	for f in tr_{1,3,5,13,15,17,25,27,29,37,39,41}*
  	do
  	    echo -e "food_sc\tblack" > dict_color
  	    echo -e "food_fat\tyellow" >> dict_color
  	    pergola_rules.py -i \${f} -m ${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color'
  	done

  	for f in tr_{7,9,11,19,21,23,31,33,35,43,45,47,48}*
  	do
  	    echo -e "food_sc\torange" > dict_color
  	    echo -e "food_fat\tblue" >> dict_color
  	    pergola_rules.py -i \${f} -m ${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color'
  	done

  	for f in tr_{1,3,5,13,15,17,25,27,29,37,39,41}*
  	do
  	    echo -e "food_sc\tblack" > dict_color
  	    echo -e "food_fat\torange" >> dict_color
  	    pergola_rules.py -i \${f} -m ${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color'
  	done
  	"""
}

process convert_bedGraph {

    publishDir = [path: {params.output}, mode: 'copy', overwrite: 'true']

  	input:
  	file ('batch_bg') from mice_files_bedGraph
  	file mapping_file_bG
  	val max from max_time.first()

  	output:   	
  	file 'tr*food*.bedGraph' into bedGraph_out
  	file 'tr*{water,sac}*.bedGraph' into bedGraph_out_drink
  	//file 'bedg_to_viz' into bedg_to_gviz
  	
  	"""  
  	# pergola_rules.py -i ${batch_bg} -m  ${mapping_file_bG} -max ${max} -f bedGraph -w 3600 -nt -e
  	pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt -e
  	#mkdir bedg_to_viz
  	#cp *.bedGraph bedg_to_viz
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
  	"""
}


exp_phases_bed_to_wr.subscribe {
    it.copyTo( "files/exp_phases.bed" )
}

bed_out.into { bed_out_gviz; bed_out_sushi; bed_out_igv }
bedGraph_out.into { bedGraph_out_gviz; bedGraph_out_sushi }

def wt = [1,3,5,13,15,17,25,27,29,37,39,41]
def wt_nic = [7,9,11,19,21,23,31,33,35,43,45,47,48]
def cb1 = [6,8,10,18,20,22,30,32,34,42,44,46]
def cb1_nic = [2,4,12,14,16,24,26,28,36,38,40]


def map_id_group = [ "wt" : [1,3,5,13,15,17,25,27,29,37,39,41],
                     "wt_nic" : [7,9,11,19,21,23,31,33,35,43,45,47,48],
                     "cb1" : [6,8,10,18,20,22,30,32,34,42,44,46],
                     "cb1_nic" : [2,4,12,14,16,24,26,28,36,38,40] ]

////////////////////////////
bed_out_igv.into { bed_out_wt; bed_out_wt_nic; bed_out_cb1; bed_out_cb1_nic; bed_out_wt_fat; bed_out_wt_nic_fat; bed_out_cb1_fat; bed_out_cb1_nic_fat }

result_dir_wt_food_sc = file("$params.output/1_wt_food_sc")
result_dir_wt_nic_food_sc = file("$params.output/2_wt_nic_food_sc")
result_dir_cb1_food_sc = file("$params.output/3_cb1_food_sc")
result_dir_cb1_nic_food_sc = file("$params.output/4_cb1_nic_food_sc")
result_dir_wt_food_fat = file("$params.output/5_wt_food_fat")
result_dir_wt_nic_food_fat = file("$params.output/6_wt_nic_food_fat")
result_dir_cb1_food_fat = file("$params.output/7_cb1_food_fat")
result_dir_cb1_nic_food_fat = file("$params.output/8_cb1_nic_food_fat")

result_dir_wt_food_sc.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_wt_food_sc"
}

result_dir_wt_food_fat.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_wt_food_fat"
}

result_dir_wt_nic_food_sc.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_wt_nic_food_sc"
}

result_dir_wt_nic_food_fat.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_wt_nic_food_fat"
}

result_dir_cb1_food_sc.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_cb1_food_sc"
}

result_dir_cb1_food_fat.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_cb1_food_fat"
}

result_dir_cb1_nic_food_sc.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_cb1_nic_food_sc"
}

result_dir_cb1_nic_food_fat.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_cb1_nic_food_fat"
}

bed_out_wt.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    wt.contains(id.toInteger()) && food == "sc"
}.subscribe {
    it.copyTo( result_dir_wt_food_sc.resolve ( it.name ) )
}

bed_out_wt_nic.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    wt_nic.contains(id.toInteger()) && food == "sc"
}.subscribe {
    it.copyTo( result_dir_wt_nic_food_sc.resolve ( it.name ) )
}

bed_out_cb1.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    cb1.contains(id.toInteger()) && food == "sc"
}.subscribe {
    it.copyTo( result_dir_cb1_food_sc.resolve ( it.name ) )
}

bed_out_cb1_nic.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    cb1_nic.contains(id.toInteger()) && food == "sc"
}.subscribe {
    it.copyTo( result_dir_cb1_nic_food_sc.resolve ( it.name ) )
}

bed_out_wt_fat.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    wt.contains(id.toInteger()) && food == "fat"
}.subscribe {
    it.copyTo( result_dir_wt_food_fat.resolve ( it.name ) )
}

bed_out_wt_nic_fat.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    wt_nic.contains(id.toInteger()) && food == "fat"
}.subscribe {
    it.copyTo( result_dir_wt_nic_food_fat.resolve ( it.name ) )
}

bed_out_cb1_fat.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    cb1.contains(id.toInteger()) && food == "fat"
}.subscribe {
    it.copyTo( result_dir_cb1_food_fat.resolve ( it.name ) )
}

bed_out_cb1_nic_fat.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    cb1_nic.contains(id.toInteger()) && food == "fat"
}.subscribe {
    it.copyTo( result_dir_cb1_nic_food_fat.resolve ( it.name ) )
}


//def wt_food_fat = [1,3,5,13,15,17,25,27,29,37,39,41]
//def wt_nic_food_fat = [1,3,5,13,15,17,25,27,29,37,39,41]
//def cb1_food_fat = [1,3,5,13,15,17,25,27,29,37,39,41]
//def cb1_nic_food_fat = [1,3,5,13,15,17,25,27,29,37,39,41]

/*
food_type = [ "sc", "fat" ]
group = ["wt", "wt_nic", "cb1", "cb1_nic"]

group.each {
    def group = it
    food_type.each {
        def food = it
        println ( "===" + group + "_" + it )
        def name_dir = "${params.output}${group}_${food}"
        //println "====== dir inside each: $name_dir"

        dir = file(name_dir)

        dir.with {
            if( !empty() ) { deleteDir() }
            mkdirs()
            println "Created: $dir"
        }
    }
}

bed_out_igv.into { bed_out_wt; bed_out_wt_nic; bed_out_cb1; bed_out_cb1_nic; bed_out_wt_fat; bed_out_wt_nic_fat; bed_out_cb1_fat; bed_out_cb1_nic_fat }

map_channels =   [ "wt_sc" : bed_out_wt,
                   "wt_nic_sc" : bed_out_wt_nic,
                   "cb1_sc" : bed_out_cb1,
                   "cb1_nic_sc" : bed_out_cb1_nic,
                   "wt_fat" : bed_out_wt_fat,
                   "wt_nic_fat" : bed_out_wt_nic_fat,
                   "cb1_fat" : bed_out_cb1_fat,
                   "cb1_nic_fat" : bed_out_cb1_nic_fat ]


for (group_mice in group) {
    println ("group=================== $group_mice")
    for (food_i in food_type) {

        println ("*************** $food_i")

        //bed_out_igv.into { bed_out_igv; bed_out_write }

        //bed_out_write.flatten().println()
        key = group_mice + "_" + food_i
        map_channels.get(key).flatten().filter {
        //bed_out_write.flatten().filter {

            def id = it.name.split("\\_")[1]

            def food = it.name.split("\\_")[4].split("\\.")[0]
            println ("group mice &&&&&& =================== $group_mice")

            map_id_group.get(group_mice).contains(id.toInteger()) && food == food_i

            //group.contains(id.toInteger()) && food == food_type

        }.subscribe {
            println ("***********" + it.name)
            it.copyTo( "${params.output}${group}_${food_i}".resolve ( it.name ) )
        }

    }
}
*/
/*
bed_out_igv.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (it.name)
    wt.contains(id.toInteger()) && food == "sc"
}
*/
/*
.subscribe {
    it.copyTo( result_dir_wt_food_sc.resolve ( it.name ) )
}
*/
/*
result_dir_wt_food_sc = file("$params.output/wt_food_sc")
result_dir_wt_food_fat = file("$params.output/wt_food_fat")
result_dir_wt_nic_food_sc = file("$params.output/wt_nic_food_sc")
result_dir_wt_nic_food_fat = file("$params.output/wt_nic_food_fat")
result_dir_cb1_food_sc = file("$params.output/cb1_food_sc")
result_dir_cb1_food_fat = file("$params.output/cb1_food_fat")
result_dir_cb1_nic_food_sc = file("$params.output/cb1_nic_food_sc")
result_dir_cb1_nic_food_fat = file("$params.output/cb1_nic_food_fat")

result_dir_wt_food_sc.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
    println "Created: $result_dir_wt_food_sc"
}

bed_out_igv.flatten().filter {
    def id = it.name.split("\\_")[1]
    def food = it.name.split("\\_")[4].split("\\.")[0]
    //println (food)
    wt.contains(id.toInteger()) && food == "sc"
}.subscribe {
    it.copyTo( result_dir_wt_food_sc.resolve ( it.name ) )
}


.subscribe {
	  it.copyTo( result_dir_heatmap.resolve ( it.name ) )
}


bed_out_igv.collect {
    def id = it.name.split("\\.")[1]

    println id
}

motion_files_flat_p = motion_files.map { name_mat, motion_f ->
                          motion_f.collect {
        	                def motion = it.name.split("\\.")[1]
                          [ it, name_mat, it.name, motion ]
                          }
                      }
                      .flatMap()
*/


//().println()

/*
motion_f.collect {
        	                def motion = it.name.split("\\.")[1]
                          [ it, name_mat, it.name, motion ]
                          }
*/

//.filter { it[4] == tag_str1 }

process gviz_visualization {

    publishDir = [path: "plot", mode: 'copy', overwrite: 'true']

    input:

    file 'exp_info' from exp_info
    file 'bed_dir/*' from bed_out_gviz.collect()
    file 'bedgr_dir/*' from bedGraph_out_gviz.collect()
    file exp_phases_bed from exp_phases_bed_gviz

    output:
    file '*.tiff' into gviz

  	"""
    mice_gviz_visualization.R --f_experiment_info=${exp_info} --path_bed_files=bed_dir --path_to_bedGraph_files=bedgr_dir --path_to_phases_file=${exp_phases_bed}
  	"""
}

process sushi_visualization {

    publishDir = [path: "plot", mode: 'copy', overwrite: 'true']

    input:

    file 'exp_info' from exp_info
    file 'bed_dir/*' from bed_out_sushi.collect()
    file 'bedgr_dir/*' from bedGraph_out_sushi.collect()
    file exp_phases_bed from exp_phases_bed_sushi

    output:
    file '*.pdf' into sushi

  	"""
    mice_sushi_visualization.R --f_experiment_info=${exp_info} --path_bed_files=bed_dir --path_to_bedGraph_files=bedgr_dir --path_to_phases_file=${exp_phases_bed}
  	"""

}
