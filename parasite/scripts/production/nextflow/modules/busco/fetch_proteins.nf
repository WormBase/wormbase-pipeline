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
  val(dbs)
  
  storeDir "${params.output_dir}/${dbs.species}/protein/"

  output:
  path "translations.fa", emit: protein
  val "${dbs.species}", emit:output_dir
  val dbs.database, emit:db_name
  val dbs.augustus_species
  val dbs.odb

  script:
  """
  mkdir -p ${params.output_dir}/${dbs.species}/protein/
  perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname ${dbs.database} -user ${params.user} -dnadbhost ${params.host} -dnadbport ${params.port} -dnadbname ${dbs.database} -dnadbuser ${params.user} -file translations.fa  ${params.dump_params}
  """
}
