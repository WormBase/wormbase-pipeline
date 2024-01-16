import re
import argparse
from ProductionMysql import *
from ProductionUtils import *

def parse_busco_output(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

        match = re.search(r'\s+C:(\d+\.\d+)%\[S:(\d+\.\d+%),D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)', content)
        if match:
            complete = float(match.group(1))
            duplicated = float(match.group(3))
            fragmented = float(match.group(4))
            missing = float(match.group(5))
            tnumber = int(match.group(6))
        else:
            exit_with_error(f"No match found in the BUSCO output file {file_path}. Exiting.")
    
    # Close file connection
    file.close()

    return complete, duplicated, fragmented, missing, tnumber

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse BUSCO output file")
    parser.add_argument("-s", "--file_path", type=str, help="Path to the BUSCO short summary file", required=True)
    parser.add_argument("-t", "--host", type=str, help="Host name of the database", required=True)
    parser.add_argument("-u", "--user", type=str, help="Username of the database", required=False)
    parser.add_argument("-r", "--port", type=str, help="Port of the database", required=False)
    parser.add_argument("-p", "--password", type=str, help="Password of the database", required=False)
    parser.add_argument("-m", "--mode", type=str, help="BUSCO mode (assembly or annotation)", required=True)
    parser.add_argument("-v", "--busco_version", type=str, help="BUSCO version", required=True)
    parser.add_argument("-d", "--database", type=str, help="Database name", required=True)
    args = parser.parse_args()

    if args.file_path:
        complete, duplicated, fragmented, missing, tnumber = parse_busco_output(args.file_path)
    else:
        print("Please provide the path to the BUSCO output file using the -s/--file_path argument.")
    
    if args.mode not in ["assembly", "annotation"]:
        print("Please provide a valid BUSCO mode using the -m/--mode argument.")
        exit(1)
    
    
    core = Core(args.host, args.database, writable=True)

    core.remove_busco_scores(mode=args.mode, busco_version=args.busco_version)

    core.add_busco_scores(mode=args.mode, busco_version=args.busco_version,
    complete=complete, duplicated=duplicated, fragmented=fragmented,
    missing=missing, tnumber=tnumber)
    