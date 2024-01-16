import os
import sys
import time
import gzip
import datetime
import glob
import subprocess
import re
import urllib.request
import requests
import paramiko
import tempfile
import shutil
import random
import validators
import sqlalchemy
import string
import json
import csv
import pandas as pd

PARASITE_FTP_URL="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases"
PARASITE_HTTP_FTP_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases"

def ftp_suffix_dict():
    return { "gfasta" : "genomic.fa.gz",
             "gfasta_masked" : "genomic_masked.fa.gz",
             "gfasta_softmasked" : "genomic_softmasked.fa.gz",
             "pfasta" : "protein.fa.gz"}


# Logging
def dtnow():
    ppdtnow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return (str(ppdtnow))


def exit_with_error(error_message):
    print(dtnow() + ": ERROR - " + error_message + ".", file=sys.stderr)
    sys.exit(1)


def check_file_exists(file, to_raise=False):
    if not os.path.isfile(file):
        if to_raise:
            exit_with_error(file + "doesn't exist or cannot be accessed.")
            raise IOError
        return (False)
    return (True)


def check_dir_exists(file, to_raise=False):
    if not os.path.isdir(file):
        if to_raise:
            exit_with_error(file + "doesn't exist or cannot be accessed.")
            raise IOError
        return (False)
    return (True)


def print_info(info_message):
    print(dtnow() + ": INFO - " + info_message + ".")


def print_warning(info_message):
    print(dtnow() + ": WARNING - " + info_message + ".")

def print_error(info_message):
    print(dtnow() + ": ERROR - " + info_message + ".")

def print_w_indent(message):
    print("\t" + message)


def pnl():
    print("\n")


def is_valid_directory_name(directory_name):
    """
    Returns True if the directory_name can be safely used as a directory name,
    otherwise returns False.
    """
    # Check if directory_name is empty or contains invalid characters
    if not directory_name or not all(char.isalnum() or char in "-_ " for char in directory_name):
        return False

    # Check if directory with same name already exists
    if os.path.exists(directory_name):
        return False

    return True


def fix_directory_name(directory_name):
    """
    Tries to fix the given directory_name so that it can be used as a directory name.
    """
    # Remove invalid characters from directory_name
    directory_name = ''.join(char for char in directory_name if char.isalnum() or char in "-_ ")

    # Remove leading/trailing spaces from directory_name
    directory_name = directory_name.strip()

    # Replace spaces and other non-alphanumeric characters with hyphens
    directory_name = '-'.join(directory_name.split())

    if not is_valid_directory_name(directory_name):
        exit_with_error("Cannot fix: "+directory_name)

    return directory_name


def create_soft_link_for_all_dir_contents(src_dir, dst_dir):
    """
    Creates a soft link to all the contents of src_dir in dst_dir.
    """
    # Get list of all files/directories in src_dir
    src_contents = os.listdir(src_dir)

    # Create soft links to all files/directories in src_dir in dst_dir
    for content in src_contents:
        src_path = os.path.join(src_dir, content)
        dst_path = os.path.join(dst_dir, content)
        os.symlink(src_path, dst_path)


def check_if_file_exists(xfile):
    if os.path.exists(xfile):
        pass
    else:
        exit_with_error("File: " + xfile + " does not exist. Exiting.")


# File Handling/Parsing
def decompress(infile, tofile):
    with open(infile, 'rb') as inf, open(tofile, 'w', encoding='utf8') as tof:
        decom_str = gzip.decompress(inf.read()).decode('utf-8')
        tof.write(decom_str)


def csvlines2list(csv_path, delimiter=",", skiplines=0):
    """Function that takes a path for a csv file, it opens it and returns each line as a
    python list element."""
    csv_in = csv_path.strip()
    with open(csv_in, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimiter)
        outlist = []
        skiplinescounter = 0
        while skiplinescounter < skiplines:
            next(csv_reader)
            skiplinescounter += 1
        for row in csv_reader:
            row = [x.strip() for x in row]
            row = [x.split("\xef\xbb\xbf")[1] if x.startswith("\xef\xbb\xbf") else x for x in row]
            outlist.append(row)
    return (outlist)


# Misc 
def getnum(text):
    retnnum = int(text.split("/")[-1].split("_")[5])
    return (retnnum)


def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


def flatten(t):
    return [item for sublist in t for item in sublist]


def getnum(text):
    retnnum = int(text.split("/")[-1].split("_")[5])
    return (retnnum)


def regex_match_dbs(pattern, databases):
    r = re.compile(".*" + pattern + ".*")
    filtdb_list = list(filter(r.match, databases))
    return (filtdb_list)


def regex_match_one_db(pattern, databases):
    r = re.compile(".*" + pattern + ".*")
    filtdb_list = list(filter(r.match, databases))
    if len(filtdb_list) == 0:
        exit_with_error("no db for: " + pattern)
    elif len(filtdb_list) > 1:
        exit_with_error("multiple_dbs_for: " + pattern)
    else:
        return (filtdb_list[0])


def url_file_exists(path):
    r = requests.head(path, stream=True)
    return r.status_code == requests.codes.ok

def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def download_if_url_or_copy_file_command(path, output=None):
    if validators.url(path):
        command = "curl -H 'Connection: keep-alive' --keepalive-time 2 " + \
                  ("--output " + output + " " if output else "") + \
                  path + ";"
    else:
        command = "rsync -avz " + \
                  path + \
                  (" " + output if output else "") + \
                  ";"
    return command

def argparse_core_db(cdb_input, staging_server, previous_staging_server):
    if cdb_input in staging_server.core_databases:
        return "{0}:{1]".format(cdb_input, staging_server.host)
    elif cdb_input in previous_staging_server.core_databases:
        return "{0}:{1}".format(cdb_input, previous_staging_server.host)
    else:
        raise ValueError

def argparse_fasta_file(fasta_input):
    if not (fasta_input.endswith(".fa")) or (fasta_input.endswith(".fasta")) \
            or (fasta_input.endswith(".fa.gz")) or (fasta_input.endswith(".fasta.gz")):
        raise ValueError
    if not ":" in fasta_input:
        raise ValueError
    return fasta_input

def parasite_release_data_genomes(parasite_data_dir, staging_server):
    return [x for x in os.listdir(parasite_data_dir) if x in staging_server.release_genomes]

def parasite_genome2core(genome, staging_server):
    core_dbs = staging_server.core_dbs(genome)
    if len(core_dbs)>1:
        exit_with_error("The are more than one core dbs in {0} for {1}".format(staging_server.host,genome))
    return core_dbs[0]

def parasite_core2genome(core_db):
    return "_".join(core_db.split("_")[0:3])

def release_from_dbname(dbname, dbtype="core"):
    """
    Extracts the release information from a given database name.

    The function expects the database name to follow a specific format with a prefix `dbtype`, followed by an underscore,
    a number with at least two digits representing the release, another underscore, and a number with at least three digits.
    It extracts the release information from the given `dbname` and returns it.

    Args:
        dbname (str): The name of the database containing the release information.
        dbtype (str, optional): The specific type of the database to look for (default: "core").

    Returns:
        str or None: The extracted release information if found, or None if not found.

    Example:
        dbname = "core_18_108"
        release = release_from_dbname(dbname)
        print(release)
        # Output: 18
    """    
    pattern = rf"{dbtype}_(\d{{2,}})_(\d{{3,}})"
    match = re.search(pattern, dbname)
    if match:
        return match.group(1)
    else:
        return None    

def parasite_core2previouscore(core_db, previous_staging_server):
    genome = parasite_core2genome(core_db)
    return parasite_genome2core(genome, previous_staging_server)

def parasite_release_data_cores(parasite_data_dir, staging_server):
    parasite_data_genomes = parasite_release_data_genomes(parasite_data_dir, staging_server)
    return [parasite_genome2core(x, staging_server) for x in parasite_data_genomes]

def get_ftp_species_url_for_version(version, parasite_ftp_url):
    ps_version = str(version)
    url = f"{parasite_ftp_url}/WBPS{ps_version}/species/"
    return url

def get_ftp_species_for_version(version, parasite_ftp_url=PARASITE_HTTP_FTP_URL):
    url = get_ftp_species_url_for_version(version, parasite_ftp_url)
    page = urllib.request.urlopen(url)
    content = page.read().decode('utf-8')
    regex = r"href=\"(\w+_\w+)/"
    genomes = re.findall(regex, content)
    return genomes

def get_previous_release_genomes_list(version, parasite_ftp_url=PARASITE_HTTP_FTP_URL):
    genomes = get_ftp_species_for_version(version, parasite_ftp_url)
    url = get_ftp_species_url_for_version(version, parasite_ftp_url)

    genomes_list = []
    for genome in genomes:
        genome_name = re.sub(r"_", " ", genome)
        bioprojects_url = f"{url}/{genome}/"
        bioprojects_page = urllib.request.urlopen(bioprojects_url)
        bioprojects_content = bioprojects_page.read().decode('utf-8')
        bioprojects_regex = r"href=\"(\w+)/"
        bioprojects = re.findall(bioprojects_regex, bioprojects_content)
        for bioproject in bioprojects:
            bioproject_name = bioproject.strip().lower()
            genome_bioproject_name = f"{genome}_{bioproject_name}"
            genomes_list.append(genome_bioproject_name)

    return genomes_list

def get_genome_ftp_file_url_from_ftp_release(genome, version, filetype, parasite_ftp_url=PARASITE_HTTP_FTP_URL):
    ps_version = str(version)
    genus_species = "_".join(genome.split("_")[0:2])
    bioproject = "_".join(genome.split("_")[2:]).upper()

    url = f"{parasite_ftp_url}/WBPS{ps_version}/species/{genus_species}/{bioproject}"

    filesuffix = ftp_suffix_dict()[filetype]

    file_url = f"{url}/{genus_species}.{bioproject}.WBPS{ps_version}.{filesuffix}"

    if url_file_exists(file_url):
        return file_url
    else:
        exit_with_error(f"{file_url} doesn't exit. Exiting.")

def get_all_genomes_ftp_file_url_from_ftp_release(version, filetype, parasite_ftp_url=PARASITE_HTTP_FTP_URL):
    release_genomes = get_previous_release_genomes_list(version)
    file_urls = []
    for genome in release_genomes:
        file_urls.append(get_genome_ftp_file_url_from_ftp_release(genome, version, filetype, parasite_ftp_url))

    return file_urls

def calculate_n50_command(fasta, outfile, assembly_stats_path="assembly-stats"):
    return f'{assembly_stats_path} {fasta} | grep "N50 =" > {outfile};'

def open_remote_file(remote_address):
    # Parse the remote address
    username, hostname, filepath = remote_address.split('@')[0], remote_address.split('@')[1].split(':')[0], remote_address.split(':')[1]

    # Create an SSH client
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()

    # Connect to the remote host
    ssh.connect(hostname=hostname, username=username)

    # Open the remote file
    sftp = ssh.open_sftp()
    remote_file = sftp.open(filepath, 'r')

    # Define a function to close the remote file
    def close_remote_file():
        remote_file.close()
        sftp.close()
        ssh.close()

    # Register the function to be called when the object is garbage collected
    remote_file.__del__ = close_remote_file

    return remote_file


def check_remote_file_exists(remote_address):
    # Parse the remote address
    username, hostname, filepath = remote_address.split('@')[0], remote_address.split('@')[1].split(':')[0], remote_address.split(':')[1]

    # Create an SSH client
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()

    try:
        # Connect to the remote host
        ssh.connect(hostname=hostname, username=username)

        # Check if the file exists
        sftp = ssh.open_sftp()
        try:
            sftp.stat(filepath)
            return True
        except FileNotFoundError:
            return False
        finally:
            sftp.close()
    finally:
        ssh.close()

def rename_remote_file(source_address, destination_address):
    # Parse the source and destination addresses
    source_username, source_hostname, source_filepath = source_address.split('@')[0], source_address.split('@')[1].split(':')[0], source_address.split(':')[1]
    destination_username, destination_hostname, destination_filepath = destination_address.split('@')[0], destination_address.split('@')[1].split(':')[0], destination_address.split(':')[1]

    # Create an SSH client
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()

    try:
        # Connect to the remote hosts
        ssh.connect(hostname=source_hostname, username=source_username)
        ssh.connect(hostname=destination_hostname, username=destination_username)

        # Rename the remote file
        sftp = ssh.open_sftp()
        try:
            sftp.rename(source_filepath, destination_filepath)
        finally:
            sftp.close()
    finally:
        ssh.close()







