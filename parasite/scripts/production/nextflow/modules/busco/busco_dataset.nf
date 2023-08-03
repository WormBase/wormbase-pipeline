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
  val odb_version

  output:
  tuple val(dbs), stdout

  script:
  """
  #!python
  from ProductionMysql import *
  import json
  import os
  core = Core('${server.host}', '${dbs.database}')

  #json_path = 'busco_datasets.${dbs.database}.json'

  #busco_dict = {
  #  'odb': core.busco_dataset('${odb_version}'),
  #  'augustus_species': core.busco_augustus_species()
  #}

  #json_object = json.dumps(busco_dict, indent=4)

  #with open(json_path, "w") as outfile:
  #  outfile.write(json_object)

  odb = core.busco_dataset('${odb_version}')
  augustus_species = core.busco_augustus_species()
  print(f"{odb};{augustus_species}")
  """

}