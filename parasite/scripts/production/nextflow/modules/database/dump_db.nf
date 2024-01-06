// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

process DUMP_DB {
    publishDir "$out_dir/coredb/$db.division", mode: 'copy'
    tag "$db.species"
    label "variable_2_8_32"
    maxForks 10

    input:
        val server
        val db
        val out_dir

    output:
        path "*.sql.gz"

    script:
        """
        db_pass=""
        if [ "${server.password}" != "" ]; then
            db_pass="--password '${server.password}'"
        fi

        mysqldump '${db.database}' \
            --host '${server.host}' \
            --port '${server.port}' \
            --user '${server.user}' \
            \$db_pass \
            | gzip > ${db.species}.sql.gz
        """
}
