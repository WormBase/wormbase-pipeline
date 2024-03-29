import sqlalchemy
from sqlalchemy import insert
from ProductionUtils import *
import subprocess
import os
import re
from ensembl.database import DBConnection, UnitTestDB

# Usage examples:
# Core(staging.host, "mansoni").ftp_filename()
# Core(previous_staging.host, "mansoni").ftp_filename()
# staging.core_databases


# Example 1: Loading this script and execute a command in a core database:
# from ProductionMysql import * 
# my_core = Core(staging.host, "schistosoma_mansoni_prjea36577_core_18_108_1")
# execution = my_core.connect().execute("SELECT * FROM gene LIMIT 10;")
# print([x for x in execution])


# Example 2: Loading this script and get a meta_value from the meta table
# of a core database:
# from ProductionMysql import *
# my_core = Core(staging.host, "schistosoma_mansoni_prjea36577_core_18_108_1")
# my_core.meta_value("species.taxonomy_id")

STAGING_HOST = os.environ['PARASITE_STAGING_MYSQL']
PRODUCTION_HOST = 'mysql-ps-prod-1'
PREVIOUS_STAGING_HOST = os.environ['PREVIOUS_PARASITE_STAGING_MYSQL']
PARASITE_VERSION = os.environ['PARASITE_VERSION']
ENSEMBL_VERSION = os.environ['ENSEMBL_VERSION']
WORMBASE_VERSION = os.environ['WORMBASE_VERSION']
WORMBASE_ENSEMBL_SCRIPTS = os.path.join(os.environ['WORM_CODE'], 'scripts', 'ENSEMBL', 'scripts')
DUMP_GENOME_SCRIPT = os.path.join(WORMBASE_ENSEMBL_SCRIPTS, "dump_genome.pl")
DUMP_PROTEIN_SCRIPT = os.path.join(WORMBASE_ENSEMBL_SCRIPTS, "dump_proteins.pl")

class Staging:
    def __init__(self, STAGING_HOST, writeable=False, ps_release=PARASITE_VERSION, e_release=ENSEMBL_VERSION):
        self.host = STAGING_HOST
        if writeable:
            if self.host.startswith("mysql-ps-"):
                if self.host.endswith(".ebi.ac.uk"):
                    split_host = self.host.split(".ebi.ac.uk")
                    self.host = split_host[0] + "-w"
                else:
                    self.host = STAGING_HOST + "-w"
            else:
                exit_with_error(f"It is not allowed to get a writeable version of the non ps sever {STAGING_HOST}.")
        
        self.url = subprocess.check_output([self.host, 'details', 'url']).strip().decode()
        self.script = subprocess.check_output([self.host, 'details', 'script']).strip().decode()
        self.engine = sqlalchemy.create_engine(self.url)
        self.insp = sqlalchemy.inspect(self.engine)
        self.db_list = self.insp.get_schema_names()
        self.core_databases = [x for x in self.db_list if "_core_" in x]
        self.variation_databases = [x for x in self.db_list if "_variation_" in x]
        self.relsuffix = "_"+ps_release+"_"+e_release
        self.release_core_databases = [x for x in self.core_databases if "_"+ps_release+"_"+e_release in x]
        self.genomes = ["_".join(x.split("_")[0:3]) for x in self.core_databases]
        self.release_genomes = ["_".join(x.split("_")[0:3]) for x in self.release_core_databases]
        self.species = list(set(["_".join(x.split("_")[0:2]) for x in self.core_databases]))
        self.release_species = list(set(["_".join(x.split("_")[0:2]) for x in self.release_core_databases]))

    def core_dbs(self, pattern):
        return (regex_match_dbs(pattern, self.core_databases))

    def variation_dbs(self, pattern):
        return (regex_match_dbs(pattern, self.variation_databases))

    def species(self, pattern):
        fdbs = self.dbs(pattern)
        return (['_'.join(x.split("_")[0:3]) for x in fdbs])

    def is_core(self, input):
        if input in self.core_dbs("_core_" + os.environ["PARASITE_VERSION"] + "_"):
            return True
        else:
            return False

    def is_a_core(self, input):
        if input in self.core_databases:
            return True
        else:
            return False

    def is_genome(self, input):
        if input in self.release_genomes:
            return True
        else:
            return False

    def is_a_genome(self, input):
        if input in self.genomes:
            return True
        else:
            return False

    def species_phylum(self, pattern):
        if "diphyllobothrium_latum" in pattern:
            pattern = "dibothriocephalus_latus_prjeb1206"
        try:
            coredb = regex_match_dbs(pattern, self.core_databases)[0]
        except IndexError:
            exit_with_error(f"Couldn't find a core db for {pattern}")

        return Core(self.host, coredb).phylum()
    
    def core_dbs_with_busco(self, version, pattern=None, release=False, opposite=False):
        if (pattern==None) and (release==False):
            all_cores = self.core_databases
        elif pattern and release==False:
            all_cores = self.core_dbs(pattern)
        elif release==True and pattern==None:
            all_cores = self.release_core_databases
        elif release==True and pattern:
            all_cores = regex_match_dbs(pattern, self.release_core_databases)
        
        core_dbs_with_buscos = []

        for dbname in all_cores:
            if Core(self.host, dbname).has_busco_score(version):
                core_dbs_with_buscos.append(dbname)
        
        return core_dbs_with_buscos


class Production:
    def __init__(self, PRODUCTION_HOST, pattern="", writeable=False):
        self.host = PRODUCTION_HOST
        if writeable:
            self.host = PRODUCTION_HOST + "-w"
        self.url = subprocess.check_output([self.host, 'details', 'url']).strip().decode()
        self.script = subprocess.check_output([self.host, 'details', 'script']).strip().decode()
        self.engine = sqlalchemy.create_engine(self.url)
        self.insp = sqlalchemy.inspect(self.engine)
        self.db_list = self.insp.get_schema_names()
        self.pattern = pattern

    def db(self):
        return (regex_match_one_db(self.pattern, self.db_list))

    def connect(self):
        dbc = DBConnection(self.url + self.db())
        return (dbc)

class Core(Staging):
    def __init__(self, STAGING_HOST, pattern, writable=False):
        super().__init__(STAGING_HOST, writable)
        self.pattern = pattern
        self.databases = self.core_databases
        self.core_url = self.url + self.db()
        self.engine = sqlalchemy.create_engine(self.core_url)
    def db(self):
        return (regex_match_one_db(self.pattern, self.databases))
    def species_name(self):
        spe_cies_bp = '_'.join(self.db().split("_")[0:3])
        return (spe_cies_bp)
    def release(self):
        return release_from_dbname(self.db(),"core")
    def connect(self):
        dbc = DBConnection(self.url + self.db())
        return (dbc)
    def meta_value(self, pattern):
        META_VALUE_SQL = f"SELECT meta_key, meta_value FROM meta WHERE meta_key LIKE '%{pattern}%';"
        result = self.connect().execute(META_VALUE_SQL).fetchall()
        if result:
            if len(result) > 1:
                return [{x[0]: x[1]} for x in result]
            else:
                return result[0][1]
        else:
            return []
    def phylum(self):
        species_classifications_dict = self.meta_value('species.classification')
        species_classifications = flatten([list(x.values()) for x in species_classifications_dict])
        if "Nematoda" in species_classifications:
            return "Nematoda"
        elif "Platyhelminthes" in species_classifications:
            return "Platyhelminthes"
        else:
            exit_with_error("Could not infer phylum from: "+" ".join(species_classifications))
    def assembly_default(self):
        return self.meta_value('assembly.default')
    def ftp_filename_n_filename_without_version(self):
        ftp_id = self.meta_value('species.ftp_genome_id')
        ftp_species = '_'.join(self.species_name().split("_")[0:2])
        return (ftp_species + "/" + ftp_id + "/" + ftp_species + "." + ftp_id + ".")
    def ftp_filename_n_filename(self):
        ftp_id = self.meta_value('species.ftp_genome_id')
        ftp_species = '_'.join(self.species_name().split("_")[0:2])
        return (ftp_species + "/" + ftp_id + "/" + ftp_species + "." + ftp_id + "." + "WBPS" + PARASITE_VERSION)
    def ftp_filename(self):
        ftp_id = self.meta_value('species.ftp_genome_id')
        ftp_species = '_'.join(self.species_name().split("_")[0:2])
        return (ftp_species + "." + ftp_id + "." + "WBPS" + PARASITE_VERSION)
    def dump_genome_command(self, outfile=False, mask=False, softmask=False, ebi_header_prefix=False):
        command = "perl " + \
                  DUMP_GENOME_SCRIPT + " " + \
                  self.script + " " + \
                  "--dbname " + self.db() + \
                  (" --outfile " + outfile if outfile else "") + \
                  (" --mask" if mask else "") + \
                  (" --softmask" if softmask else "") + \
                  ";"
        return command
    def dump_protein_command(self, outfile=False, canonical_only=True):
        command = "perl " + \
                  DUMP_PROTEIN_SCRIPT + " " + \
                  self.script + " " + \
                  "--dbname " + self.db() + \
                  (" --outfile " + outfile if outfile else "") + \
                  (" --canonical_only" if canonical_only else "") + \
                  ";"
        return command
    def has_busco_score(self, version):
        busco_scores = self.meta_value(f"busco{version}")
        if len(busco_scores)>=10:
            return True
        else:
            return False
    def busco_dataset(self, odb_version):
        if self.phylum()=="Nematoda":
            return f"nematoda_odb{odb_version}"
        elif self.phylum()=="Platyhelminthes":
            return f"metazoa_odb{odb_version}"
    def busco_augustus_species(self):
        if self.phylum()=="Nematoda":
            return f"caenorhabditis"
        elif self.phylum()=="Platyhelminthes":
            return f"schistosoma"
    def remove_busco_scores(self, mode, busco_version):
        META_VALUE_SQL = f"DELETE FROM meta WHERE meta_key LIKE '{mode}.busco{busco_version}%';"
        self.connect().execute(META_VALUE_SQL)
    def add_busco_scores(self, mode, busco_version, complete, duplicated, fragmented, missing, tnumber):
        busco_dict = {f"{mode}.busco{busco_version}_complete":f"{complete}",
                      f"{mode}.busco{busco_version}_duplicated":f"{duplicated}",
                      f"{mode}.busco{busco_version}_fragmented":f"{fragmented}",
                      f"{mode}.busco{busco_version}_missing":f"{missing}",
                      f"{mode}.busco{busco_version}_number":f"{tnumber}"}
        for key, value in busco_dict.items():
            META_VALUE_SQL = f"INSERT INTO meta (meta_key, meta_value) VALUES ('{key}', '{value}');"
            self.connect().execute(META_VALUE_SQL)
    def remove_omark_scores(self):
        META_VALUE_SQL = f"DELETE FROM meta WHERE meta_key LIKE 'omark.%';"
        self.connect().execute(META_VALUE_SQL)
    def add_omark_scores(self, single, duplicated, missing, accurate, partial, fragmented, contamination, unknown):
        omark_dict = {"omark.single":f"{single}",
                      "omark.duplicated":f"{duplicated}",
                      "omark.missing":f"{missing}",
                      "omark.consistent":f"{accurate}",
                      "omark.partial":f"{partial}",
                      "omark.fragmented":f"{fragmented}",
                      "omark.contamination":f"{contamination}",
                      "omark.unknown":f"{unknown}"}
        for key, value in omark_dict.items():
            META_VALUE_SQL = f"INSERT INTO meta (meta_key, meta_value) VALUES ('{key}', '{value}');"
            self.connect().execute(META_VALUE_SQL)

class Variation(Staging):
    def __init__(self, STAGING_HOST, pattern, writable=False):
        super().__init__(STAGING_HOST, writable)
        self.pattern = pattern
        self.databases = self.variation_databases

    def db(self):
        return (regex_match_one_db(self.pattern, self.databases))

    def species_name(self):
        spe_cies_bp = '_'.join(self.db().split("_")[0:3])
        return (spe_cies_bp)

    def release(self):
        return release_from_dbname(self.db(),"variation")
    
    def core_dbname(self):
        return (Core(STAGING_HOST, self.species_name()).db())

    def connect(self):
        dbc = DBConnection(self.url + self.db())
        return (dbc)

    def meta_value(self, pattern):
        META_VALUE_SQL = "SELECT meta_key, meta_value FROM meta WHERE meta_key LIKE \"%{0}%\";".format(pattern)
        get_meta_value_res = [x for x in self.connect().execute(META_VALUE_SQL)]
        if len(get_meta_value_res) > 1:
            return ({x[0]: x[1] for x in get_meta_value_res})
        else:
            return ([x[1] for x in get_meta_value_res][0])

    def ftp_filename_n_filename(self):
        ftp_id = self.meta_value('species.ftp_genome_id')
        ftp_species = '_'.join(self.species_name().split("_")[0:2])
        return (ftp_species + "/" + ftp_id + "/" + ftp_species + "." + ftp_id + "." + "WBPS" + PARASITE_VERSION)

    def ftp_filename(self):
        ftp_id = self.meta_value('species.ftp_genome_id')
        ftp_species = '_'.join(self.species_name().split("_")[0:2])
        return (ftp_species + "." + ftp_id + "." + "WBPS" + PARASITE_VERSION)

staging = Staging(STAGING_HOST)
previous_staging = Staging(PREVIOUS_STAGING_HOST)
stagingw = Staging(STAGING_HOST, writeable=True)
previous_stagingw = Staging(STAGING_HOST, writeable=True)

def is_parasite_genome(pattern, staging=staging):
    if len(staging.core_dbs(pattern)) >= 1:
        return True
    else:
        return False

def core_which_staging(core_db, staging=staging, previous_staging=previous_staging):
    if staging.is_a_core(core_db):
        this_staging=staging
    elif previous_staging.is_a_core(core_db):
        this_staging=previous_staging
    else:
        exit_with_error("Sorry but " + core_db + " is not a core db in staging or previous staging.")
    return(this_staging)
