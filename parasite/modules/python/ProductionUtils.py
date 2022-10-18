import os
import sys
import time
import gzip
import datetime
import glob
import subprocess
import re
import requests
import random
import sqlalchemy
import string
import csv
import pandas as pd


# Logging
def dtnow():
    ppdtnow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return (str(ppdtnow))


def exit_with_error(error_message):
    print(dtnow() + ": ERROR - " + error_message + ".")
    raise ValueError


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


def print_w_indent(message):
    print("\t" + message)


def pnl():
    print("\n")


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
    r = requests.head(path)
    return r.status_code == requests.codes.ok


def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


