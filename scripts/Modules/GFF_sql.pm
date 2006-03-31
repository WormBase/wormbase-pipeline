#!/usr/bin/env perl

###########################
# class to manage GFF stuff
# by Michael Han
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-03-31 09:25:04 $
######

package GFF_sql;

#use warnings;
use strict;
use DBI;

sub new {
    my $class   = shift;
    my ($options) = @_;
    my $self    = {};
    bless $self, $class;
    $self->{debug} = 1 if $options->{'-debug'};
    $self->{build} = 1 if $options->{'-build'};
    $self->{fallback} =1 if $options->{'-fallback'};
    print "debug => ", $self->{'debug'}, "\n" if $options->{'-debug'};
    $self->initialize;
    return $self;
}

sub initialize {
    my $self = shift;
    my ( $load, $chromosome ) = @_;
    
    my $dsn="DBI:mysql:database=mh6;host=mcs2a;port=3316";
    my $user='mh6';
    my $pass='mh6';
    if ($self->{'build'}) { $dsn='DBI:mysql:database=mh6_build;host=mcs2a;port=3316' }
    elsif ($self->{fallback}) { $dsn='DBI:mysql:database=worm_testdb;host=ecs1f';$user='wormadmin';$pass='worms'}
 
#    $self->{dbh} = DBI->connect( "DBI:mysql:database=worm_testdb;host=ecs1f",
#        'wormadmin', 'worms' )
     $self->{dbh} = DBI->connect( $dsn,$user,$pass) || die "Cannot connect: $DBI::errstr";
    return 1;
}

########################
# clean(chromosome_name)
#
sub clean {
	my $self=shift;
	my $chromosome=shift;
	$self->{dbh}->do("DROP TABLE $chromosome");
	$self->{dbh}->do("CREATE TABLE $chromosome (tag_id SMALLINT UNSIGNED NOT NULL, type_id SMALLINT UNSIGNED NOT NULL, start INT UNSIGNED NOT NULL, stop INT UNSIGNED NOT NULL, frame ENUM('.','0','1','2') NOT NULL default '.', orientation ENUM('+','-','.') NOT NULL default '.',fluff TEXT not null, INDEX (start) )");

}


######################
# load_gff(file_name,chromosome_name)
# load GFF into db
#
sub load_gff {

    use IO::File;

    my $self = shift;
    my ( $file, $chromosome, $overwrite ) = @_;
    $self->{dbh}->{AutoCommit} = 0;

    #prepare db
    if ($overwrite) {
    
        $self->{dbh}->do("DROP TABLE IF EXISTS $chromosome");
        $self->{dbh}->do(
"CREATE TABLE $chromosome (tag_id SMALLINT UNSIGNED NOT NULL, type_id SMALLINT UNSIGNED NOT NULL, start INT UNSIGNED NOT NULL, stop INT UNSIGNED NOT NULL, frame ENUM('.','0','1','2') NOT NULL default '.', orientation ENUM('+','-','.') NOT NULL default '.',fluff TEXT not null, INDEX (start) )"
        );
   }
    my $sth =
      $self->{dbh}->prepare(
"INSERT INTO $chromosome (tag_id,type_id,start,stop,frame,orientation,fluff) VALUES (?,?,?,?,?,?,?)"
      );
    $self->{dbh}->commit;
   my $tagh = $self->{dbh}->prepare("SELECT id  FROM gff_tag WHERE name=?");
   my $typeh= $self->{dbh}->prepare("SELECT id  FROM gff_types WHERE name=?");
    #loop over GFF and split it
    my $fh = new IO::File $file;
    die "cannot open $file\n" if not( defined $fh );
    while (<$fh>) {
        my @fields      = split;
        my $chromosome_ = shift @fields;          #never used
        my $tag         = shift @fields;
        my $type        = shift @fields;
        my $start       = shift @fields;
        my $stop        = shift @fields;
        my $score       = shift @fields;
        my $orientation = shift @fields;
        my $frame       = shift @fields;
        my $note        = join( ' ', @fields );
        next if $_ =~ /^\#/;
	$typeh->  execute($type);
	$tagh ->  execute($tag);
	my ($tag_sql) = $tagh ->fetchrow_array;
        my ($type_sql)= $typeh->fetchrow_array;
        $sth->execute( $tag_sql, $type_sql, $start, $stop, $frame, $orientation,
            $note );
    }
    $fh->close;
    $self->{dbh}->commit;
}

##########
# generate_tags
# generate tag and types from a gff

sub generate_tags {
    use IO::File;

    my $self = shift;
    my ($file) = @_;
    my %tags;
    my %types;

    $self->{dbh}->{AutoCommit} = 0;

    #prepare db

    my $typeh =
      $self->{dbh}->prepare("INSERT IGNORE INTO gff_types (name) VALUES (?)");
    my $tagh = $self->{dbh}->prepare("INSERT IGNORE INTO gff_tag (name) VALUES (?)");
    $self->{dbh}->commit;

    #loop over GFF and split it
    my $fh = new IO::File $file;
    die "cannot open $file\n" if not( defined $fh );
    while (<$fh>) {
        my @fields      = split;
        my $chromosome_ = shift @fields;          #never used
        my $tag         = shift @fields;
        my $type        = shift @fields;
        my $start       = shift @fields;
        my $stop        = shift @fields;
        my $frame       = shift @fields;
        my $orientation = shift @fields;
        my $no_clue     = shift @fields;
        my $note        = join( ' ', @fields );
        next if $_ =~ /^\#/;
        $tags{$tag}   = 0;
        $types{$type} = 0;
    }
    foreach my $key ( keys %types ) { $typeh->execute($key) }
    foreach my $key ( keys %tags )  { $tagh->execute($key) }

    $fh->close;
    $self->{dbh}->commit;
}

##########
# get_gff(chromosome,start,stop,@limits)
#
sub get_gff {
    my $self = shift;
    my ( $chromosome, $start, $stop, @limit ) = @_;
    my $sth;
    my $queryprefix =
"SELECT gff_tag.name as feature ,gff_types.name as source,start,stop,frame,orientation,fluff FROM $chromosome LEFT JOIN gff_types ON $chromosome.type_id=gff_types.id LEFT JOIN gff_tag ON $chromosome.tag_id=gff_tag.id WHERE ";

    if (@limit) {
        $sth =
          $self->{dbh}->prepare( "$queryprefix gff_types.name IN (\'"
              . ( join( "\',\'", @limit ) )
              . "\') AND((start BETWEEN ? AND  ?) OR (stop BETWEEN ? AND  ?) OR (start < ? AND stop > ?))"
          )
          || die "cannot prepare statement:  $DBI::errstr";
    }
    else {

        $sth =
          $self->{dbh}->prepare(
"$queryprefix (start BETWEEN ? AND  ?) OR (stop BETWEEN ? AND  ?) OR (start < ? AND stop > ?)"
          )
          || die "cannot prepare statement:  $DBI::errstr";

    }

    $sth->execute( $start, $stop, $start, $stop, $start, $stop );
    my $tbl_ary_ref = $sth->fetchall_arrayref( {} );
    return @$tbl_ary_ref;
}

##########
# get_gff2(chromosome,start,stop,%limit)
#
sub get_gff2 {
    my $self = shift;
    my ( $chromosome, $start, $stop, $limit ) = @_;
    $limit->{'start'} = $start;
    $limit->{'stop'}  = $stop;
    return $self->get_chr( $chromosome, $limit );
}

#######
# get_by_source($chromosome,$type)
# 2bremoved
sub get_by_source {
    my $self = shift;
    my ( $chromosome, $type ) = @_;
    return $self->get_chr( $chromosome, { 'source' => $type } );
}

######
# get_bestblat_byfluff($chromosome,$fluff)
# hardcoded to BLAT_mRNA_BEST (debug)
#
sub get_bestblat_byfluff {
    my $self = shift;
    my ( $chromosome, $fluff ) = @_;

    return $self->get_chr( $chromosome,
        { 'feature' => 'BLAT_mRNA_BEST', 'fluff' => $fluff } );
}

##########
# get_chr($chr,$id)
#
sub get_chr {
    my $self = shift;
    my ( $chromosome, $limit ) = @_;

    my @where;
    my $queryprefix =
"SELECT gff_tag.name as feature ,gff_types.name as source ,start,stop,frame,orientation,fluff FROM $chromosome LEFT JOIN gff_types ON $chromosome.type_id=gff_types.id LEFT JOIN gff_tag ON $chromosome.tag_id=gff_tag.id";
    $queryprefix = "$queryprefix WHERE " if $limit;
    if ( $limit && $limit->{'start'} && $limit->{'stop'}) {
        push @where, " (stop >= $limit->{'start'} AND  start <= $limit->{'stop'}) ";
    }

    if ( $limit->{'feature'} ) {
        push @where, " gff_tag.name=\'" . $limit->{'feature'} . "\' ";
    }
    if ( $limit->{'source'} ) {
        push @where, " gff_types.name=\'" . $limit->{'source'} . "\' ";
    }
    if ( $limit->{'fluff'} ) {
        push @where, " fluff LIKE \'\%$limit->{'fluff'}\%\' ";
    }

    my $flat_where = join( 'AND', @where );
    print "$queryprefix $flat_where\n" if $self->{'debug'};
    my $sth = $self->{dbh}->prepare("$queryprefix $flat_where") || die "cannot prepare statement:  $DBI::errstr";
    #print STDERR "...... $queryprefix $flat_where\n";
    $sth->execute();
    my $tbl_ary_ref = $sth->fetchall_arrayref( {} );
    map { $_->{'chromosome'} = $chromosome } @$tbl_ary_ref;
    return @$tbl_ary_ref;
}

#############
# getall
sub getall{
    my $self = shift;
    my ($limits) = @_;
    my @all;
    my @chr = $self->get_chr_tables();
    foreach my $chromosome (@chr) {
        my @hits = $self->get_chr( $chromosome, $limits );
        push @all, @hits if scalar(@hits) >= 1;
    }
    return @all;
}


#############
# getall_oftype
#
sub getall_oftype {
    my $self = shift;
    my ($type) = @_;
    my @all;

    my @chr = $self->get_chr_tables();
    foreach my $chromosome (@chr) {
        my @hits = $self->get_chr( $chromosome, { 'source' => $type } );
        push @all, @hits if scalar(@hits) >= 1;
    }
    return @all;
}

############
# get_allele_byid($id)
#

sub get_allele_byid {
    my $self = shift;
    my ($id) = @_;
    my @all;

    my @chr = $self->get_chr_tables();
    foreach my $chromosome (@chr) {
        my @hits =
          $self->get_chr( $chromosome,
            { 'fluff' => $id, 'feature' => 'Allele' } );
        map { $_->{chromosome} = $chromosome } @hits;
        return $hits[0] if scalar(@hits) == 1;
    }
    return undef;
}

#######
# get_best_blat($id)
# only cDNAs at the moment
sub get_best_blat {
    my $self = shift;
    my ($id) = @_;

    my @chr = $self->get_chr_tables();
    foreach my $chromosome (@chr) {
        my @hits = $self->get_bestblat_byfluff( $chromosome, $id );
        map { $_->{chromosome} = $chromosome } @hits;
        return @hits if scalar(@hits) >= 1;
    }
    return undef;
}

#######
# get chromosome tables as array
#
sub get_chr_tables {
    my $self = shift;
    my @chrs;
    my $sth = $self->{dbh}->table_info( '', '', '', 'TABLE' );
    my $tables = $sth->fetchall_arrayref;
    foreach my $tbl (@$tables) {
        my $flat_tb = "@{$tbl}";    # needs fixing for undef array elements
        if ( $flat_tb =~ /(CHROMOSOME_\w+)\sTABLE/ ) {
            push @chrs, $1;
        }
    }
    return @chrs;
}

############
# make sure the connection is closed upon destruction
#
sub DESTROY {
    my $self = shift;
    $self->{dbh}->disconnect if $self->{dbh};
}

1;
__END__

=head1 NAME

GFF_sql..pm - Class to access GFF in mySQL

=head1 SYNOPSIS

use GFF_sql;

to connect to the database:

my $object=Gff_sql->new();

=head3 to load a GFF into the database:

$object->load_gff('mygff.file','CHROMOSOME_I',$overwrite);

add $overwrite to parameters to drop table before inserting/updating data.

=head3 to get all exons and introns of chromosome 1 between 100 and 200 bp

my @hits=$object->get_gff('CHROMOSOME_I',100,200,'exons','introns');

=head3 to get all genes on chromosome 1:

my @genes=$object->get_by_source('CHROMOSOME_I','gene')

=head3 to sort hits:

@sorted_hits = sort { $a->{'start'} <=> $b->{'start'}} @unsorted_hits

=head1 DESCRIPTION

Class to access and generate GFF annotatino to the chromosomes

=head3 Database:

    chromosome_name(start,stop,orientation,fluff) <- one per chromosome

    gff_types      (sources)

    gff_tag        (features)

    partially normalized to increase speed, without any indices yet.


=head3 load_gff($filename,$chromosome_name)

returns: 1 on success else it should die

needs an GFF filename and the internal name of the chromosome
loads the GFF into the database after dropping the old table.


=head3 get_gff($chromosome_name,$start,$stop,@limits)

returns an array of all lines. Lines will be returned if they overlapp with tha area specified by start-stop. In addition@limits restricts the search to certain tags.
The type returnes is an array of hashrefs with the column names as keys (source,feature.name,start,stop,frame,orientation,fluff).

=head3 get_by_source($chromosome,$type)

returns an array of all lines. Lines will be returened if they contain the tag.
The type returned is an array of hashrefs with the column names as keys (feature,source,start,stop,frame,orientation,fluff)

=head3 get_chr($chromosome,{'source'=>x ,'feature'=>y,'fluff'=>z, 'start' => a, 'stop' => b})

This is the method most used and most flexible.

gets all hits from a chromosome base on the limits set by the hash.
returns an arry of hashrefs.
The "fluff" is used as "LIKE %$fluff%" SQLwise.

=head3 get_all({'source'=>x ,'feature'=>y,'fluff'=>z, 'start' => a, 'stop' => b})

like get_chr, just iterates over all chromosomes.

=head3 get_chr_tables()

returns an array of the names of the tables containing the GFF information per chromosome.

=head3 getall_oftype($tupe)

returns an array of all entries of a given type found in the database.

=head3 get_allele_byid($id)

returns lines from the db as %

=head3 DESTROY

Destructor to make sure that the database connection is closed on death of the object.


=head3 Hash

the GFF lines form the database are stored in a hash which is normally stored in an array and returned from methods.

%line={

    'start'          => start coordinate on the chrmomosome 
    'stop'           => stop coordinate
    'orientation'    => orientation as +/-/.
    'fluff'          => description from GFF
    'source'         => type like  intron/exon/gene/...
    'feature'        => tag like  curated/... 
    'frame'          => frame 0/1/2 
    'chromosome'     => chromosome (like CHROMOSOME_I)

}

=head2 oneliners:

to load all autoace GFFs:

including new tags and overwriting of the old table
perl -mGFF_sql -e 'my $a=GFF_sql->new; foreach $i (@ARGV){$i=~/(CHROMOSOME_\w+)\.gff/;my $name=$1;$a->generate_tags($i);$a->load_gff($i,$name,1)}' /wormsrv2/autoace/CHROMOSOMES/*.gff

show all chromosome tables:

perl -mGFF_sql -e 'my $a=GFF_sql->new();my @line=$a->get_chr_tables;foreach $i(@line){print "$i\t"}'

get all genes from CHROMOSOME_I and sort them by stop:

perl -mGFF_sql -e '$bla=GFF_sql->new();@line=$bla->get_by_source($ARGV[0],"gene");@sorted = sort { $a->{"stop"} <=> $b->{"stop"} } @line;foreach my $i(@sorted){ foreach my $key(keys %$i){print $i->{$key},"\t"}print "\n"}' CHROMOSOME_I

=head1 AUTHOR

$Author: mh6 $

=head1 VERSION

$Date: 2006-03-31 09:25:04 $

=cut
