#!/usr/bin/env nextflow
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// // dealing with custom $NXF_WORK based workDir
// //   don't do anything if "-work-dir (-w)" option specified on command line
// include { updateWorkDirAsNeeded } from '../../modules/utils/utils.nf'
// updateWorkDirAsNeeded("dumper_pipeline2")

nextflow.enable.dsl=2

// default params
params.help = false
params.prefix = ''
params.dbname_re = ''
params.output_dir = './dumper_output'
params.password = ''

// Print usage
def helpMessage() {
  log.info '''
        Mandatory arguments:
        --host, --port, --user         Connection parameters to the SQL servers we getting core db(s) from

        Optional arguments:
        --password                     Password part of the connection parameters
        --prefix                       Core dabase(s) name prefixes
        --dbname_re                    Regexp to match core db name(s) against
        --output_dir                   Name of Output directory to gather prepared outfiles. Default -> 'Output_GenomePrepare'.
        --mode                         Busco mode: genome or protein, default is to run both.
        --help                         This usage statement.
        '''
}

// Check mandatory parameters
if (params.help) {
    helpMessage()
    exit 0
}

def create_server(params) {
    server = [
        "host": params.host,
        "port": params.port,
        "user": params.user,
        "password": ""
    ]
    if (params.password) {
        server["password"] = params.password
    }
    return server
}

def create_filter_map(params) {
    filter_map = [
        "prefix": "",
        "dbname_re": ""
    ]
    if (params.prefix) {
        filter_map["prefix"] = params.prefix
    }
    if (params.dbname_re) {
        filter_map["dbname_re"] = params.dbname_re
    }
    return filter_map
}

if (params.host && params.port && params.user && params.output_dir) {
    server = create_server(params)
    filter_map = create_filter_map(params)
} else {
    exit 1, "Missing server parameters"
}

busco_mode = []
if (params.mode instanceof java.lang.String) {
  busco_mode = [params.mode]
}
else {
  busco_mode = params.mode
}

if (params.output_dir) {
    outdir = params.output_dir
}
else {
    exit 1, "Missing output directory"
}

// include { DUMP_SQL } from '../../subworkflows/dump_sql/main.nf'
// include { DUMP_METADATA } from '../../subworkflows/dump_metadata/main.nf'
include { BUSCO_ODB } from '../../modules/busco/busco_dataset.nf'
include { BUSCO_AUGUSTUS_SPECIES } from '../../modules/busco/busco_dataset.nf'
include { DB_FACTORY } from '../../modules/database/db_factory.nf'
include { read_json } from '../../modules/utils/utils.nf'
include { FETCH_GENOME } from '../../modules/busco/fetch_genome.nf'
include { BUSCO_GENOME_RUN } from '../../modules/busco/busco_genome_run.nf'

// Run main workflow
workflow {
    
    dbs = DB_FACTORY(server, filter_map)
        .map(it -> read_json(it))
        .flatten()
    
    buscoModes = Channel.fromList(busco_mode)

    dbs_w_odb = BUSCO_ODB(dbs, server, params.odb_version) 
                | map { [database: it[0].database, species: it[0].species, division: it[0].division, odb: it[1]]}

    busco_dbs = BUSCO_AUGUSTUS_SPECIES(dbs_w_odb, server) 
                | map { [database: it[0].database, species: it[0].species, division: it[0].division, odb: it[0].odb, augustus_species: it[1]] }
    
    if (busco_mode.contains('genome')) {
        genome_data = FETCH_GENOME(busco_dbs)
        BUSCO_GENOME_RUN(genome_data)
    }

    // DUMP_SQL(server, dbs, filter_map, params.output_dir)
    // DUMP_METADATA(server, dbs, filter_map, params.output_dir)
}
