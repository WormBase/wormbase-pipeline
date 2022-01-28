#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python

import sys
import os
import re

args = sys.argv

core_db = args[1]

if len(args)<2:
    print("Usage: python get_bioproject.py core_db\nExample:  python get_bioproject.py ancylostoma_ceylanicum_prjna231479_core_17_105_1")
    raise ValueError

try:
    bioproject = re.search("(prj[edbna]+\d+)_core", core_db).group(1)
    print(bioproject)
except AttributeError:
    print("")


