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
PREVIOUS_STAGING_HOST = os.environ['PREVIOUS_PARASITE_STAGING_MYSQL']
PARASITE_VERSION = os.environ['PARASITE_VERSION']
ENSEMBL_VERSION = os.environ['ENSEMBL_VERSION']
WORMBASE_VERSION = os.environ['WORMBASE_VERSION']
WORMBASE_ENSEMBL_SCRIPTS = os.path.join(os.environ['WORM_CODE'], 'scripts', 'ENSEMBL', 'scripts')
DUMP_GENOME_SCRIPT = os.path.join(WORMBASE_ENSEMBL_SCRIPTS, "dump_genome.pl")


class Staging:
    def __init__(self, STAGING_HOST, writeable=False):
        self.host = STAGING_HOST
        if writeable:
            self.host = STAGING_HOST + "-w"
        self.url = subprocess.check_output([self.host, 'details', 'url']).strip().decode()
        self.script = subprocess.check_output([self.host, 'details', 'script']).strip().decode()
        self.engine = sqlalchemy.create_engine(self.url)
        self.insp = sqlalchemy.inspect(self.engine)
        self.db_list = self.insp.get_schema_names()
        self.core_databases = [x for x in self.db_list if "_core_" in x]
        self.variation_databases = [x for x in self.db_list if "_variation_" in x]

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

class Core(Staging):
    def __init__(self, STAGING_HOST, pattern, writable=False):
        super().__init__(STAGING_HOST, writable)
        self.pattern = pattern
        self.databases = self.core_databases

    def db(self):
        return (regex_match_one_db(self.pattern, self.databases))

    def species(self):
        spe_cies_bp = '_'.join(self.db().split("_")[0:3])
        return (spe_cies_bp)

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
        ftp_species = '_'.join(self.species().split("_")[0:2])
        return (ftp_species + "/" + ftp_id + "/" + ftp_species + "." + ftp_id + "." + "WBPS" + PARASITE_VERSION)

    def ftp_filename(self):
        ftp_id = self.meta_value('species.ftp_genome_id')
        ftp_species = '_'.join(self.species().split("_")[0:2])
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



class Variation(Staging):
    def __init__(self, STAGING_HOST, pattern, writable=False):
        super().__init__(STAGING_HOST, writable)
        self.pattern = pattern
        self.databases = self.variation_databases

    def db(self):
        return (regex_match_one_db(self.pattern, self.databases))

    def species(self):
        spe_cies_bp = '_'.join(self.db().split("_")[0:3])
        return (spe_cies_bp)

    def core(self):
        return (Core(STAGING_HOST, self.species()).db())

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
        ftp_species = '_'.join(self.species().split("_")[0:2])
        return (ftp_species + "/" + ftp_id + "/" + ftp_species + "." + ftp_id + "." + "WBPS" + PARASITE_VERSION)

    def ftp_filename(self):
        ftp_id = self.meta_value('species.ftp_genome_id')
        ftp_species = '_'.join(self.species().split("_")[0:2])
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
