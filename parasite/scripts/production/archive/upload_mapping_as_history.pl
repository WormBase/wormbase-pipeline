#!/usr/bin/env perl
# Copy the tables from previous db to start with a clean slate. Then
# upload an ID mapping, and add kill and creation events for everything else.
# We only add events for genes, while Ensembl also tracks individual exons and transcripts.
use strict;
use warnings;
use ProductionMysql;
use SpeciesFtp;
use IO::Uncompress::Gunzip;
use Digest::MD5 qw/md5_hex/;
use Getopt::Long;
use ProductionMysql;
use MysqlPair;
my $db_command = "$ENV{PARASITE_STAGING_MYSQL}-w";
my $previous_db_command = "$ENV{PREVIOUS_PARASITE_STAGING_MYSQL}";
my %args;
GetOptions (
  'species=s'             => \$args{species}, # can be a pattern
  'previous_species=s'    => \$args{previous_species}, # can be a pattern
  'db_command=s'          => \$db_command,
  'previous_db_command=s' => \$previous_db_command,
  'mapping=s'             => \$args{mapping},
  'proteins=s'            => \$args{proteins},
);
die "Required --mapping <mapping file>" unless -f $args{mapping};
die "Required --species <species_pattern>" unless $args{species};
my $mysql = ProductionMysql->new($db_command);
my $previous_mysql = ProductionMysql->new($previous_db_command);
$args{core_db} = $mysql->core_db($args{species});
$args{previous_core_db} = $previous_mysql->core_db( $args{previous_species} );
my @core_db_words = split "_", $args{core_db};
$args{parasite_version} = $core_db_words[4];
$args{assembly_name} = $mysql->meta_value($args{core_db}, "assembly.name");

my @previous_core_db_words = split "_", $args{previous_core_db};
$args{previous_parasite_version} = $previous_core_db_words[4];
$args{previous_assembly_name} = $previous_mysql->meta_value($args{previous_core_db}, "assembly.name");

$args{proteins} //= SpeciesFtp->release($args{previous_parasite_version})->path_to($previous_mysql->species($args{previous_species}), "protein.fa");
die "Usage: --proteins ftp_path.fa.gz. Not a file: $args{proteins}" unless -s $args{proteins};

my $dbc = $mysql->dbc($args{core_db});
my $transcript_adaptor = $mysql->adaptor($args{core_db}, 'Transcript');

my $mysql_pair = MysqlPair->new($previous_mysql, $mysql);
# Restore mappings. This is usually a clean wipe, unless processing an update of an update.
for my $table (qw/mapping_session stable_id_event gene_archive peptide_archive/){
  $mysql_pair->dump($args{previous_species}, $args{species}, $table);
}

my $mapping_session_id = &start_mapping_session($dbc, %args);
&archive_where_seq_changed($dbc, $transcript_adaptor, $mapping_session_id, %args);

# populate stable_id_event
&add_link_events_for_mapped_genes($dbc, $mapping_session_id, %args);
&add_link_events_for_unmapped_genes_matching_ids($dbc, $mapping_session_id, %args);
&add_kill_events_for_archived_and_not_current_or_mapped_genes($dbc, $mapping_session_id, %args);
&add_creation_events_for_current_and_not_archived_or_mapped_genes($dbc, $mapping_session_id, %args);

sub start_mapping_session {
    my ($dbc, %args) = @_;
    my $insert_mapping_session_dbh = $dbc->prepare("insert into mapping_session(old_db_name,new_db_name,old_release,new_release,old_assembly,new_assembly,created) values(?,?,?,?,?,?,NOW());");
    $insert_mapping_session_dbh->execute(@args{qw/previous_core_db core_db previous_parasite_version parasite_version previous_assembly_name assembly_name/});
    my $get_session_id_dbh = $dbc->prepare("select max(mapping_session_id) from mapping_session");
    $get_session_id_dbh->execute;
    my ($mapping_session_id, @others) = $get_session_id_dbh->fetchrow_array;
    die unless $mapping_session_id;
    die if @others;
    $insert_mapping_session_dbh->finish;
    $get_session_id_dbh->finish;
    return $mapping_session_id;
}
sub archive_where_seq_changed {
    my ($dbc, $transcript_adaptor, $mapping_session_id, %args) = @_;
    my $insert_1_dbh = $dbc->prepare("insert into peptide_archive (md5_checksum, peptide_seq) values (?,?);");
    my $insert_2_dbh = $dbc->prepare("insert into gene_archive (gene_stable_id, gene_version, transcript_stable_id, translation_stable_id, peptide_archive_id, mapping_session_id) values (?, 1, ?, ?, (select max(peptide_archive_id) from peptide_archive), $mapping_session_id);");
    local $/='>';
    my $z = IO::Uncompress::Gunzip->new($args{proteins});
    while(<$z>){
        chomp;
        my ($l, @seq) = split "\n";
        next unless @seq;
        my $seq = join "", @seq;
        my ($translation, $transcript, $gene) = $l =~ /(.*) transcript(?:_id)?=(.*) gene(?:_id)?=(.*)$/;
        my $transcript_currently = $transcript_adaptor->fetch_by_stable_id($transcript);
        my $translation_currently = $transcript_currently ? $transcript_currently->translate : undef;
        my $seq_currently = $translation_currently ? $translation_currently->seq : undef;
        next if defined $seq_currently and $seq_currently eq $seq;
        $insert_1_dbh->execute(md5_hex($seq), $seq),
        $insert_2_dbh->execute($gene, $transcript, $translation);
    }
    $dbc->do('update gene_archive set transcript_version = 0, translation_version = 0;');
}
sub add_link_events_for_mapped_genes {
  my ($dbc, $mapping_session_id, %args) = @_;
  my $current_id_present_dbh = $dbc->prepare("select 1 from gene where stable_id=? limit 1");
  my $previous_id_present_dbh = $dbc->prepare("select 1 from gene_archive where gene_stable_id=? and mapping_session_id=$mapping_session_id limit 1");
  my $add_link_event_dbh = $dbc->prepare("insert into stable_id_event (old_stable_id, old_version, new_stable_id, new_version, mapping_session_id, type, score) values (?, 1 , ?, 1, $mapping_session_id,\"gene\", ?);");
  open(my $MAPPING, "<", $args{mapping}) or die $!;
  while(<$MAPPING>){
     chomp;
     my ($previous_id, $current_id, $score) = split "\t";
     $previous_id_present_dbh->execute($previous_id);
     warn "Previous ID in mapping but not in archive: $previous_id" and next unless $previous_id_present_dbh->fetchrow_array;
     $current_id_present_dbh->execute($current_id);
     warn "Current ID in mapping but not in archive: $current_id" and next unless $current_id_present_dbh->fetchrow_array;
     $score//=1;
     $add_link_event_dbh->execute($previous_id, $current_id, $score);
  }
  close $MAPPING;
}
sub add_link_events_for_unmapped_genes_matching_ids {
  my ($dbc, $mapping_session_id, %args) = @_;
  my $get_archived_ids_remaining_with_no_events_dbh = $dbc->prepare("select distinct gene_stable_id from gene_archive left join stable_id_event on (gene_archive.gene_stable_id = stable_id_event.old_stable_id ) join gene on (gene_archive.gene_stable_id = gene.stable_id) where gene_archive.mapping_session_id=$mapping_session_id and stable_id_event.mapping_session_id=$mapping_session_id and stable_id_event.old_stable_id is null");
  my $add_unknown_link_event_dbh = $dbc->prepare("insert into stable_id_event (old_stable_id, old_version, new_stable_id, new_version, mapping_session_id, type, score) values (?, 1, ?, 1, $mapping_session_id, \"gene\", 0);");
  $get_archived_ids_remaining_with_no_events_dbh->execute;
  while (my ($archived_id) = $get_archived_ids_remaining_with_no_events_dbh->fetchrow_array) {
     $add_unknown_link_event_dbh->execute($archived_id,$archived_id);
  }
}
sub add_kill_events_for_archived_and_not_current_or_mapped_genes {
  my ($dbc, $mapping_session_id, %args) = @_;
  my $get_archived_ids_not_remaining_with_no_events_dbh= $dbc->prepare("select distinct gene_stable_id from gene_archive left join stable_id_event on (gene_archive.gene_stable_id = stable_id_event.old_stable_id) left join gene on (gene_archive.gene_stable_id = gene.stable_id) where gene_archive.mapping_session_id=$mapping_session_id and stable_id_event.mapping_session_id=$mapping_session_id and stable_id_event.old_stable_id is null and gene.stable_id is null");
  my $add_kill_event_dbh = $dbc->prepare("insert into stable_id_event (old_stable_id, old_version, new_stable_id, new_version, mapping_session_id, type, score) values (?, 1, NULL, NULL, $mapping_session_id, \"gene\", 0);");
  $get_archived_ids_not_remaining_with_no_events_dbh->execute;
  while (my ($archived_id) = $get_archived_ids_not_remaining_with_no_events_dbh->fetchrow_array) {
    $add_kill_event_dbh->execute($archived_id);
  }
}
sub add_creation_events_for_current_and_not_archived_or_mapped_genes {
  my ($dbc, $mapping_session_id, %args) = @_;
  my $get_stable_ids_with_no_events_dbh = $dbc->prepare("select stable_id from gene left join stable_id_event on ( gene.stable_id = stable_id_event.new_stable_id ) where stable_id_event.new_stable_id is null");
  my $add_creation_event_dbh = $dbc->prepare("insert into stable_id_event (old_stable_id, old_version, new_stable_id, new_version, mapping_session_id, type, score) values (NULL, NULL , ?, 1, $mapping_session_id, \"gene\", 1);");
  $get_stable_ids_with_no_events_dbh->execute;
  while (my ($new_id) = $get_stable_ids_with_no_events_dbh->fetchrow_array){
    $add_creation_event_dbh->execute($new_id);
  }
}
