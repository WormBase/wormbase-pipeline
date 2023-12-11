import re
import sys
import pandas as pd

def flatten(t):
    return [item for sublist in t for item in sublist]

# Code to names mapping 
# #S:Single:S, D:Duplicated[U:Unexpected,E:Expected],M:Missing
# #A:Accurate - placements in consistent lineage[P:Partial hits,F:Fragmented], I: Inconsistent placements[P:Partial hits,F:Fragmented], C: Likely Contamination[P:Partial hits,F:Fragmented], U: Unknown
code2names = {
    'S': 'Single',
    'D': 'Duplicated',
    'M': 'Missing',
    'A': 'Accurate',
    'I': 'Inconsistent',
    'C': 'Contamination',
    'U': 'Unknown',
    'P': 'Partial',
    'F': 'Fragmented'
}

headers = ["Single", "Duplicated", "Missing", "Accurate", "Inconsistent", "Contamination", "Unknown"]

if len(sys.argv)<2:
    print("Usage: python omark_parser.py <proteins.sum file> or <header (to only print the output file headers)")
    sys.exit()
if sys.argv[1].endswith("proteins.sum") or sys.argv[1]=="header":
    pass
else:
    print("Usage: python omark_parser.py <protein.sums file> or <header (to only print the output file headers)")
    sys.exit()

if sys.argv[1] == "header":
    print("\t".join(headers))
    sys.exit()

file_path = sys.argv[1]
hog_pattern = r'^(S:\d+\.\d+?)%,(D:\d+\.\d+?)%\[\S+,\S+\],(M:\d+\.\d+?)%'
accuracy_pattern = r'^(A:\d+\.\d+?)%\[\S+,\S+\],(I:\d+\.\d+?)%\[\S+,\S+\],(C:\d+\.\d+?)%\[\S+,\S+\],(U:\d+\.\d+?)%'

with open(file_path,'r') as file:
    file_content = file.read()

hog_matches = flatten(re.findall(hog_pattern, file_content, re.MULTILINE))
acc_matches = flatten(re.findall(accuracy_pattern, file_content, re.MULTILINE))

if len(hog_matches) != 3 and len(acc_matches) != 5:
    raise Exception(f"Could not get statistics for {file_path}")

stats = hog_matches + acc_matches
stats = [x.split(":")[1] for x in stats]

print("\t".join(stats))
