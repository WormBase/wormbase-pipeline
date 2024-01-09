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
        --host, --port, --user         Required: Connection parameters to the SQL servers we getting core db(s) from

        Optional arguments:
        --password                     Required: Password part of the connection parameters
        --dbname_re                    Required: Regexp to match core db name(s) against. For multiple use the | delimiter.
        --output_dir                   Name of Output directory to gather prepared outfiles. Default -> 'Output_GenomePrepare'.
        --type                         Optional: Run mode (BUSCO/OMArk). It can be busco or omark or both, default is to run both. 
        --mode                         Optional: Busco mode. It can be genome or protein or both, default is to run both for BUSCO. For OMArk is only protein.
        --augustus_config_path         Required: Path to the Augustus/config directory (For ParaSite this is: ${PARASITE_SOFTWARE}/Augustus/config).
        --enscode                      Required: Path to the Ensembl code directory (For ParaSite this is: ${ENSEMBL_CVS_ROOT_DIR}.
        --wormpython                   Optional: Path to ParSite's pythonpath: This should be ${PYTHONPATH}.
        --wormcode                     Optional: Path to ParSite's wormcode: This should be ${WORM_CODE}.
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

run_type = []
if (params.type instanceof java.lang.String) {
  run_type = [params.type]
}
else {
  run_type = params.type
}

if (params.output_dir) {
    outdir = params.output_dir
}
else {
    exit 1, "Missing output directory"
}

if (params.busco_lineages) {
    busco_lineages = params.busco_lineages
}
else {
    exit 1, "Missing the busco lineages directory"
}

// include { DUMP_SQL } from '../../subworkflows/dump_sql/main.nf'
// include { DUMP_METADATA } from '../../subworkflows/dump_metadata/main.nf'
include { BUSCO_ODB } from '../../modules/busco/busco_dataset.nf'
include { BUSCO_AUGUSTUS_SPECIES } from '../../modules/busco/busco_dataset.nf'
include { DB_FACTORY } from '../../modules/database/db_factory.nf'
include { read_json } from '../../modules/utils/utils.nf'
include { FETCH_GENOME } from '../../modules/busco/fetch_genome.nf'
include { FETCH_PROTEINS } from '../../modules/busco/fetch_proteins.nf'
include { BUSCO3_GENOME_RUN } from '../../modules/busco/busco_genome_run.nf'
include { BUSCO3_PROTEIN_RUN } from '../../modules/busco/busco_protein_run.nf'
include { OMAMER_HOG } from '../../modules/omark/omamer_hog.nf'
include { RUN_OMARK } from '../../modules/omark/run_omark.nf'
include { BUSCOGENOME2DATABASE } from '../../modules/utils/parse_output.nf'
include { BUSCOPROTEIN2DATABASE } from '../../modules/utils/parse_output.nf'
include { OMARK2DATABASE } from '../../modules/utils/parse_output.nf'

// Run main workflow
workflow{
    
    dbs = DB_FACTORY(server, filter_map)
        .map(it -> read_json(it))
        .flatten()
    
    buscoModes = Channel.fromList(busco_mode)

    dbs_w_odb = BUSCO_ODB(dbs, server, params.odb_version) 
                | map { [database: it[0].database, species: it[0].species, division: it[0].division, odb: it[1]]}

    busco_dbs = BUSCO_AUGUSTUS_SPECIES(dbs_w_odb, server) 
                | map { [database: it[0].database, species: it[0].species, division: it[0].division, odb: it[0].odb, augustus_species: it[1]] }
    
    if (run_type.contains('busco') && busco_mode.contains('genome')) {
            genome_data = FETCH_GENOME(busco_dbs)
            bgenome_run = BUSCO3_GENOME_RUN(genome_data, busco_lineages)
            BUSCOGENOME2DATABASE(bgenome_run, server)
        }
    if (busco_mode.contains('protein') || run_type.contains('omark')) {
            protein_data = FETCH_PROTEINS(busco_dbs)
            if (run_type.contains('busco')) {
                bprotein_run = BUSCO3_PROTEIN_RUN(protein_data, busco_lineages)
                BUSCOPROTEIN2DATABASE(bprotein_run, server)
            }
            if (run_type.contains('omark')) {
                omhog = OMAMER_HOG(protein_data)
                omrun = RUN_OMARK(omhog)
                OMARK2DATABASE(omrun, server)
            }   
        }
    }
