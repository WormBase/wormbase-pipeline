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

// dump unmasked dna sequences from core db 
process FETCH_GENOME {

  label 'fetch_file'

  input:
  val(dbs)

  storeDir "${params.output_dir}/${dbs.species}/genome/"

  output:

  path "genome_toplevel.fa", emit:fasta
  val "${dbs.species}", emit:output_dir
  val dbs.database, emit:db_name
  val dbs.augustus_species
  val dbs.odb

  script:
  """
  mkdir -p ${params.output_dir}/${dbs.species}/genome/
  perl ${params.enscode}/ensembl-analysis/scripts/sequence_dump.pl -dbhost ${params.host} -dbport ${params.port} -dbname ${dbs.database} -dbuser ${params.user} -coord_system_name toplevel -toplevel -onefile -nonref -filename genome_toplevel.fa
  """

}
