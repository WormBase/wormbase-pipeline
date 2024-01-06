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

process DB_FACTORY {
    tag "DB_factory"
    label 'default'
    time '5min'

    input:
        val server
        val filter_map

    output:
        path "dbs.json"

    script:
        """
        dbname_re=''
        if [ -n '${filter_map.dbname_re}' ]; then
            dbname_re="--dbname_re ${filter_map.dbname_re}"
        fi
        db_factory --host '${server.host}' \
            --port '${server.port}' \
            --user '${server.user}' \
            --password '${server.password}' \
            --prefix '${filter_map.prefix}' \
            \$dbname_re \
            --output_json dbs.json
        """
}
