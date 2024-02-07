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
import groovy.json.JsonSlurper

def get_species_name(name) {
  final m = name.tr('[A-Z]','[a-z]').tr('.','v').split('/')[0]
  return m
}

def get_gca(name) {
  final m = name.tr('[A-Z]', '[a-z]').tr('.', 'v').replaceAll("_","").split('/').getAt(1)
  return m
}

def concatString(string1, string2, string3){
 return string1 + '_'+string2 + '_'+string3
}

def get_key_list(dict) {
    // Add quotes around each key of the dictionary to make the list compatible with Bash
    return "['" + dict.keySet().join("','") + "']"
}

def read_json(json_path) {
    slurp = new JsonSlurper()
    json_file = file(json_path)
    text = json_file.text
    // unfortunately
    //   return slurp.parseText(text)
    // doesn't work for a single element list, we suspect lazy eval
    // symptom: instead of `[a:..., b:...]` we see the same stuff in the curly brackets `{a:..., b:...}`
    not_a_lazy_val = slurp.parseText(text)
    return not_a_lazy_val
}

def updateWorkDirAsNeeded(dir_name) {
    // dealing with custom $NXF_WORK based workDir
    //   don't do anything if "-work-dir (-w)" option specified on command line
    cmd_line = binding.variables.workflow.commandLine
    cmd_line_has_wd = cmd_line.contains(" -w ") || cmd_line.contains(" -work_dir ")
    if (!cmd_line_has_wd) {
        log.info " no -work-dir (-w) option specified. Trying to build one base on NXF_WORK env"
        nxf_work_env = binding.getVariable('NXF_WORK')
        if (nxf_work_env) {
            session.workDir = ("${nxf_work_env}/dumper_pipeline" as Path).complete()
            session.workDir.mkdirs()
            workDir = session.workDir as String
        }
    }
}
