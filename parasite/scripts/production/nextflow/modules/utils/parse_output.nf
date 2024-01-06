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

// Get Busco output into the DataBase meta table
process BUSCOGENOME2DATABASE {
  scratch false

  label 'default'
   
  input:
  path summary_file
  val outdir
  val db
  val server

  script:
  """
  python ${params.wormcode}/parasite/scripts/production/core-creation/busco_parse_and_update.py -s ${summary_file} -t ${server.host} -m assembly -v ${params.busco_version_for_db} -d ${db}
  """
}

process BUSCOPROTEIN2DATABASE {
  scratch false

  label 'default'
   
  input:
  path summary_file
  val outdir
  val db
  val server

  script:
  """
  python ${params.wormcode}/parasite/scripts/production/core-creation/busco_parse_and_update.py -s ${summary_file} -t ${server.host} -m annotation -v ${params.busco_version_for_db} -d ${db}
  """
}

process OMARK2DATABASE {
  scratch false

  label 'default'
   
  input:
  path summary_file
  val outdir
  val db
  val server

  script:
  """
  python ${params.wormcode}/parasite/scripts/production/core-creation/omark_parse_and_update.py -s ${summary_file} -t ${server.host} -d ${db}
  """
}