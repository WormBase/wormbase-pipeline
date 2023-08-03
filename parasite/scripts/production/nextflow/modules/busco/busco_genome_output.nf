/*
 See the NOTICE file distributed with this work for additional information
 regarding copyright ownership.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

// Import utility functions
include { get_species_name } from './utils.nf'
include { get_gca } from './utils.nf'
include { concatString } from './utils.nf'

process BUSCO_GENOME_OUTPUT {
     
     //rename busco summary file in <production name>_gca_genome_busco_short_summary.txt
     
     label 'default'

     input:
     val outdir

     publishDir "${params.outDir}/${outdir}/",  mode: 'copy'

     script:
     """
     mkdir -p  ${params.outDir}/${outdir}/statistics
     sed  -i '/genebuild/d' ${params.outDir}/${outdir}/genome/short_summary*
     mv -f ${params.outDir}/${outdir}/genome/short_summary* ${params.outDir}/${outdir}/statistics/${concatString(get_species_name("${outdir.trim()}"),get_gca("${outdir.trim()}"),'genome_busco_short_summary.txt')}
     """
 }
