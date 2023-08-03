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

// run Busco in protein mode 
process BUSCO_PROTEIN_LINEAGE {

  label 'busco'

  input:
  file translations
  val outdir
  val db
  val busco_dataset

  output:
  path "fasta/*.txt", emit: summary_file
  val outdir, emit:species_outdir

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.outDir}/${outdir}/",  mode: 'copy'

  script:
  """
  busco -f -i ${translations}  --mode proteins -l ${busco_dataset} -c ${task.cpus} -o fasta --offline --download_path ${params.download_path}
  """
}

