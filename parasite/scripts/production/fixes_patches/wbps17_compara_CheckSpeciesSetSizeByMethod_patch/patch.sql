## Fix ORTHOLOGUES species_set_header.names ##

UPDATE species_set_header
        JOIN
    (SELECT
        *
    FROM
        (SELECT
        ssh.species_set_id,
            CONCAT(UPPER(SUBSTRING(ss1.name, 1, 1)), SUBSTRING(SUBSTRING_INDEX(SUBSTRING_INDEX(ss1.name, '_', 2), '_', - 1), 1, 3), '-', UPPER(SUBSTRING(ss2.name, 1, 1)), SUBSTRING(SUBSTRING_INDEX(SUBSTRING_INDEX(ss2.name, '_', 2), '_', - 1), 1, 3)) AS new_name
    FROM
        species_set_header ssh
    CROSS JOIN (SELECT
        *
    FROM
        species_set
    JOIN genome_db USING (genome_db_id)) ss1 USING (species_set_id)
    CROSS JOIN (SELECT
        *
    FROM
        species_set
    JOIN genome_db USING (genome_db_id)) ss2 USING (species_set_id)
    WHERE
        size = 2 AND ss1.name != ss2.name
    ORDER BY species_set_id , new_name ASC) t1
    GROUP BY species_set_id) t2 ON t2.species_set_id = species_set_header.species_set_id
SET
    name = new_name;

## Fix PARALOGUES species_set_header.names ##

UPDATE species_set_header
        JOIN
    (SELECT
        ssh.species_set_id,
            CONCAT(UPPER(SUBSTRING(gdb.name, 1, 1)), SUBSTRING(SUBSTRING_INDEX(SUBSTRING_INDEX(gdb.name, '_', 2), '_', - 1), 1, 3)) AS new_name
    FROM
        species_set_header ssh
    JOIN species_set ss USING (species_set_id)
    JOIN genome_db gdb USING (genome_db_id)
    WHERE
        size = 1) t1 ON t1.species_set_id = species_set_header.species_set_id
SET
    name = new_name;

## Fix method_link_species_set.names ##

UPDATE method_link_species_set
        JOIN
    (SELECT
    ssh.species_set_id,
        CONCAT(ssh.name, ' orthologues') AS new_name
    FROM
        method_link_species_set mlss
    JOIN species_set_header ssh USING (species_set_id)
    WHERE
        method_link_id = 201) t1 ON method_link_species_set.species_set_id = t1.species_set_id
SET
    name = new_name;

UPDATE method_link_species_set
        JOIN
    (SELECT
    ssh.species_set_id,
        CONCAT(ssh.name, ' paralogues') AS new_name
    FROM
        method_link_species_set mlss
    JOIN species_set_header ssh USING (species_set_id)
    WHERE
        method_link_id = 202) t1 ON method_link_species_set.species_set_id = t1.species_set_id
SET
    name = new_name;