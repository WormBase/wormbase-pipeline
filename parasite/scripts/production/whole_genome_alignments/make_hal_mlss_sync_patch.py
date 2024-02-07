#!/usr/bin/env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""An ad-hoc script to help synchronise Plants HAL MLSS data.

This script compares HAL MLSS data between a source and destination
database, and if there are any differences, it generates an SQL file
which may be applied to the destination database so that its HAL MLSS
data is in sync with that of the source database.

If HAL MLSS data is already in sync between the two databases,
no SQL file is generated.
"""

import argparse
from sqlalchemy import create_engine, text
import sys


major, minor, *_unused = sys.version_info
if major < 3 or (major == 3 and minor < 8):
    raise RuntimeError(f"script '{__file__}' requires Python 3.8+, exiting.")


cactus_method_types = ["CACTUS_HAL", "CACTUS_HAL_PW"]

table_to_key = {
    "genome_db": ["genome_db_id"],
    "method_link": ["method_link_id"],
    "method_link_species_set": ["method_link_species_set_id"],
    "species_set_header": ["species_set_id"],
    "species_set": ["species_set_id", "genome_db_id"],
}

rel_table_names = list(table_to_key)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source-url", required=True,
                    help="URL of synchronisation source database.")
parser.add_argument("--dest-url", required=True,
                    help="URL of synchronisation destination database.")
parser.add_argument("--release", required=True, type=int,
                    help="Current Ensembl release.")
parser.add_argument("--output-file", required=True,
                    help="Output patch file.")

args = parser.parse_args()


data = {}


data["src"] = {table_name: {} for table_name in rel_table_names}
engine = create_engine(args.source_url, future=True)
print("Checking source database ...")

query1 = text("""SELECT *
                 FROM method_link
                 WHERE type IN :cactus_method_types""")
with engine.connect() as conn:
    results = conn.execute(query1, {"cactus_method_types": cactus_method_types})
    for mapping in results.mappings():
        method_link_id = mapping["method_link_id"]
        data["src"]["method_link"][(method_link_id,)] = mapping

method_link_ids = list(data["src"]["method_link"])

query2 = text("""SELECT *
                 FROM method_link_species_set
                 WHERE method_link_id IN :method_link_ids
                 AND first_release IS NOT NULL AND first_release <= :release
                 AND (last_release IS NULL OR last_release >= :release)""")
with engine.connect() as conn:
    results = conn.execute(query2, {"method_link_ids": method_link_ids, "release": args.release})
    for mapping in results.mappings():
        mlss_id = mapping["method_link_species_set_id"]
        data["src"]["method_link_species_set"][(mlss_id,)] = mapping

mlss_ids = list(data["src"]["method_link_species_set"])

species_set_ids = [x["species_set_id"] for x in data["src"]["method_link_species_set"].values()]
species_set_ids = list(set(species_set_ids))

# We're assuming here that all the relevant species sets are
# current, based on the fact that their MLSSes are current.
query3 = text("""SELECT *
                 FROM species_set_header
                 WHERE species_set_id IN :species_set_ids""")
with engine.connect() as conn:
    results = conn.execute(query3, {"species_set_ids": species_set_ids})
    for mapping in results.mappings():
        species_set_id = mapping["species_set_id"]
        data["src"]["species_set_header"][(species_set_id,)] = mapping

query4 = text("""SELECT *
                 FROM species_set
                 WHERE species_set_id IN :species_set_ids""")
with engine.connect() as conn:
    results = conn.execute(query4, {"species_set_ids": species_set_ids})
    for mapping in results.mappings():
        species_set_id = mapping["species_set_id"]
        genome_db_id = mapping["genome_db_id"]
        data["src"]["species_set"][(species_set_id, genome_db_id)] = mapping

query5 = text("""SELECT DISTINCT genome_db.*
                 FROM genome_db
                 WHERE first_release IS NOT NULL AND first_release <= :release
                 AND (last_release IS NULL OR last_release >= :release)""")
with engine.connect() as conn:
    results = conn.execute(query5, {"release": args.release})
    for mapping in results.mappings():
        genome_db_id = mapping["genome_db_id"]
        data["src"]["genome_db"][(genome_db_id,)] = mapping


data["dst"] = {table_name: {} for table_name in rel_table_names}
engine = create_engine(args.dest_url, future=True)
print("Checking destination database ...")
dst_db_name = engine.url.database

query6 = text("""SELECT *
                 FROM method_link
                 WHERE method_link_id IN :method_link_ids""")
with engine.connect() as conn:
    results = conn.execute(query6, {"method_link_ids": method_link_ids})
    for mapping in results.mappings():
        method_link_id = mapping["method_link_id"]
        data["dst"]["method_link"][(method_link_id,)] = mapping

query7 = text("""SELECT *
                 FROM method_link_species_set
                 WHERE method_link_species_set_id IN :mlss_ids""")
with engine.connect() as conn:
    results = conn.execute(query7, {"mlss_ids": mlss_ids})
    for mapping in results.mappings():
        mlss_id = mapping["method_link_species_set_id"]
        data["dst"]["method_link_species_set"][(mlss_id,)] = mapping

with engine.connect() as conn:
    results = conn.execute(query3, {"species_set_ids": species_set_ids})
    for mapping in results.mappings():
        species_set_id = mapping["species_set_id"]
        data["dst"]["species_set_header"][(species_set_id,)] = mapping

with engine.connect() as conn:
    results = conn.execute(query4, {"species_set_ids": species_set_ids})
    for mapping in results.mappings():
        species_set_id = mapping["species_set_id"]
        genome_db_id = mapping["genome_db_id"]
        data["dst"]["species_set"][(species_set_id, genome_db_id)] = mapping

with engine.connect() as conn:
    results = conn.execute(query5, {"release": args.release})
    for mapping in results.mappings():
        genome_db_id = mapping["genome_db_id"]
        data["dst"]["genome_db"][(genome_db_id,)] = mapping


print("Comparing source and destination data ...")
inserts = []
updates = []
for table_name in rel_table_names:
    key_names = table_to_key[table_name]
    src_table = data["src"][table_name]
    dst_table = data["dst"][table_name]
    for primary_key, src_row in src_table.items():
        if primary_key in dst_table:
            dst_row = dst_table[primary_key]
            if dst_row != src_row:
                assignments = []
                for col_name, col_value in src_row.items():
                    col_value = "NULL" if col_value is None else repr(col_value)
                    assignment = f"{col_name} = {col_value}"
                    assignments.append(assignment)

                conditions = []
                for key_name, key_value in zip(key_names, primary_key):
                    assert key_value is not None
                    condition = f"{key_name} = {repr(key_value)}"
                    conditions.append(condition)

                update = f"UPDATE {table_name} SET {', '.join(assignments)} WHERE {' AND '.join(conditions)};"
                updates.append(update)
        else:
            col_values = []
            for col_value in src_row.values():
                col_value = "NULL" if col_value is None else repr(col_value)
                col_values.append(col_value)
            insert = f"INSERT INTO {table_name} VALUES ({', '.join(col_values)});"
            inserts.append(insert)


if inserts or updates:
    print(f"Writing patch to file '{args.output_file}' ...")
    with open(args.output_file, "w") as out_file_obj:
        use_stmt = f"USE {dst_db_name};"
        print(use_stmt, file=out_file_obj)

        for update in updates:
            print(update, file=out_file_obj)
        out_file_obj.write("\n")

        for insert in inserts:
            print(insert, file=out_file_obj)

    print(f"Patch written to file '{args.output_file}'.")
else:
    print("HAL MLSS data in sync between source and destination databases, no patch needed.")
