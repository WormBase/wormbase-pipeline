#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python

import sys
import os
import re

args = sys.argv

core_db = args[1]

if len(args)<2:
    print("Usage: python get_bioproject.py core_db\nExample:  python get_bioproject.py ancylostoma_ceylanicum_prjna231479_core_17_105_1")
    raise ValueError

bioproject_raw = core_db.split("_")[2]

if re.match("^prj[edbna]+\d+", bioproject_raw):
    bioproject = bioproject_raw
else:
    try:
        bioproject_match = re.search("^([^_]+)(prj[edbna]+\d+)", bioproject_raw)
        bioproject = bioproject_match.group(1) + "_" + bioproject_match.group(2)
    except AttributeError:
        print(core_db+"-----------ERROR-----------")

print(bioproject)