CREATE DATABASE IF NOT EXISTS wbps_id_minting;
USE wbps_id_minting;

CREATE TABLE organism (
    organism_id int(12) UNSIGNED AUTO_INCREMENT NOT NULL PRIMARY KEY,
    taxonomy_id int(10),
    name varchar(128),
    scientific_name varchar(128),
    is_current tinyint(1) NOT NULL DEFAULT '1',
    UNIQUE KEY name (name)
);

CREATE TABLE core_db (
    core_db_id int(12) UNSIGNED AUTO_INCREMENT NOT NULL PRIMARY KEY,
    organism_id int(10),
    parasite_version int(2),
    name varchar(128),
    is_current tinyint(1) NOT NULL DEFAULT '1',
    UNIQUE KEY name (name),
    KEY organism_index (organism_id),
    FOREIGN KEY (organism_id) REFERENCES organism(organism_id)
);

CREATE TABLE locus (
    locus_id int(12) UNSIGNED AUTO_INCREMENT NOT NULL PRIMARY KEY,
    stable_id varchar(128),
    organism_id int(12),
    is_current tinyint(1) NOT NULL DEFAULT '1',
    KEY stable_id_idx (stable_id),
    KEY organism_index (organism_id),
    FOREIGN KEY (organism_id) REFERENCES organism(organism_id)
);

CREATE TABLE gene (
    gene_id int(12) UNSIGNED AUTO_INCREMENT NOT NULL PRIMARY KEY,
    locus_id int(12),
    stable_id varchar(128),
    core_db_id int(12),
    is_current tinyint(1) NOT NULL DEFAULT '1',
    KEY stable_id_idx (stable_id),
    KEY locus_index (locus_id),
    FOREIGN KEY (core_db_id) REFERENCES core_db(core_db_id),
    FOREIGN KEY (locus_id) REFERENCES locus(locus_id)
);

CREATE TABLE transcript (
    transcript_id int(12) UNSIGNED AUTO_INCREMENT NOT NULL PRIMARY KEY,
    gene_id int(12),
    stable_id varchar(128),
    is_current tinyint(1) NOT NULL DEFAULT '1',
    KEY gene_index (gene_id),
    KEY stable_id_idx (stable_id),
    FOREIGN KEY (gene_id) REFERENCES gene(gene_id)
);




