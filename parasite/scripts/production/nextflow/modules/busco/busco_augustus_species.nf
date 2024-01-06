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

// Get Busco dataset using NCBI taxonomy in meta table 
process BUSCO_DATASET {
  scratch false

  label 'default'
   
  input:
  val dbs
  val server
  
  storeDir "${params.output_dir}/${dbs.species}/"
  
  output:
  path "busco_augustus_species.json", emit: busco_augustus_species

  script:
  """
  #!python
  from ProductionMysql import *
  import json
  core = Core('${server.host}', '${dbs.database}')
  result = core.busco_augustus_species()

  # Prepare data to be written to JSON file
  data = {
      "busco_augustus_species_result": result
  }

  # Writing data to a JSON file
  with open('busco_augustus_species.json', 'w') as json_file:
      json.dump(data, json_file, indent=4)
  """
}