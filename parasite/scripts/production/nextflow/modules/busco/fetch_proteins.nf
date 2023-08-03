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

// Dump canonical translations 
process FETCH_PROTEINS {

  label 'fetch_file'

  input:
  tuple val(species_dir),val(db), val(busco_dataset), val(mode)

  storeDir "${params.outDir}/${species_dir}/fasta/"

  output:
  path "${db}_translations.fa", emit: fasta
  val species_dir, emit: output_dir
  val db, emit:db_name
  val busco_dataset, emit:busco_dataset

  script:
  """
  perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname $db -user ${params.user} -dnadbhost ${params.host} -dnadbport ${params.port} -dnadbname $db -dnadbuser ${params.user} -file ${db}_translations.fa  ${params.dump_params}
  """
}
