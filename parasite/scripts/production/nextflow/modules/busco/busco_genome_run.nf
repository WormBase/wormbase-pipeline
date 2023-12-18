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

// run Busco in genome mode 
process BUSCO_GENOME_RUN {

  label 'busco_run'

  input:

  file genome
  val outdir
  val db
  val augustus_species
  val odb

  output:

  path "genome/*.txt", emit: summary_file
  val outdir, emit:species_outdir

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.output_dir}/${outdir}/",  mode: 'copy'

  script:
  """
  busco --lineage_dataset ${odb} --mode genome --augustus --augustus_species ${augustus_species} --cpu ${task.cpus} --offline --download_path ${params.download_path} -f -i ${genome} -o genome 
  """
}

// run Busco in genome mode 
process BUSCO3_GENOME_RUN {

  label 'busco3_run'

  input:

  file genome
  val outdir
  val db
  val augustus_species
  val odb
  val busco_lineages

  output:

  path "genome/*.txt", emit: summary_file
  val outdir, emit:species_outdir

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.output_dir}/${outdir}/",  mode: 'copy'

  script:
  """
  run_BUSCO.py -sp ${augustus_species} -l ${busco_lineages}/${odb} -o genome -i ${genome} -m genome -c ${task.cpus} -f -r
  """
}