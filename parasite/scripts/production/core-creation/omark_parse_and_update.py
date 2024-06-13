import argparse
import re

from ProductionMysql import *
from ProductionUtils import *


def parse_omark_output(file_path):
    """
    #  S:Single,
    #  D:Duplicated[U:Unexpected,E:Expected],
    #  M:Missing

    #  A:Consistent (taxonomically)[P:Partial hits,F:Fragmented],
    #  I:Inconsistent (taxonomically)[P:Partial hits,F:Fragmented],
    #  C:Likely Contamination[P:Partial hits,F:Fragmented],
    #  U:Unknown 
    """
    hog_pattern = r'^(S:\d+\.\d+?)%,(D:\d+\.\d+?)%\[\S+,\S+\],(M:\d+\.\d+?)%'
    accuracy_pattern = r'^(A:\d+\.\d+?)%\[(P:\d+\.\d+?)%,(F:\d+\.\d+?)%\],(I:\d+\.\d+?)%\[\S+,\S+\],(C:\d+\.\d+?)%\[\S+,\S+\],(U:\d+\.\d+?)%'

    with open(file_path,'r') as file:
        content = file.read()
    
    hog_matches = flatten(re.findall(hog_pattern, content, re.MULTILINE))
    acc_matches = flatten(re.findall(accuracy_pattern, content, re.MULTILINE))

    if len(hog_matches) != 3 and len(acc_matches) != 7:
        exit_with_error(f"Could not parse OMArk output file {file_path}. Exiting.")
    
    stats = hog_matches + acc_matches
    return [float(x.split(":")[1]) for x in stats]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse OMArk output file")
    parser.add_argument("-s", "--file_path", type=str, help="Path to the OMArk short summary file", required=True)
    parser.add_argument("-t", "--host", type=str, help="Host name of the database", required=True)
    parser.add_argument("-d", "--database", type=str, help="Database name", required=True)
    args = parser.parse_args()

    if args.file_path:
        single, duplicated, missing, \
        accurate, partial, fragmented, \
        inconsistent, contamination, unknown = parse_omark_output(args.file_path)
    else:
        print("Please provide the path to the OMArk output file using the -s/--file_path argument.")
    
    consistent_complete = accurate - (partial + fragmented)

    core = Core(args.host, args.database, writable=True)

    core.remove_omark_scores()

    core.add_omark_scores(single, duplicated, missing,
                          accurate, consistent_complete, partial, fragmented,
                          inconsistent, contamination, unknown)
