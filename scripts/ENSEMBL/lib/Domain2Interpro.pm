#!/usr/local/bin/perl -w
#
# Originally written by Gary Williams (gw3@sanger), modifying a script by Marc Sohrmann (ms2@sanger.ac.uk)
#
# Dumps InterPro protein motifs from ensembl mysql (protein) database to an ace file
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2015-07-09 09:39:11 $
=pod

=head1 NAME

Domain2Interpro

=head1 SYNOPSIS

use Domain2Interpro;
$m=new Domain2Interpro;
$m->get_mapping($hashref->{method}->arrayreff(array($hid,$start,$end,$hstart,$hend,$score,$evalue)))

=head1 DESCRIPTION

generalisation of Gary's Interpro dumper

=cut

package Domain2Interpro;

use strict;
####################
# Class Variables
#

=head2 class variable

       %method_database: holds the mapping of logic_name to interpro database identifiers

=cut

# define the Database names that InterPro uses in interpro.xml
# and the logic names (as specified in @methods) that search those databases
my %method_database = (
		       'scanprosite' => 'PROSITE',
		       'prints'      => 'PRINTS',
		       'pfscan'      => 'PROFILE',
		       'blastprodom' => 'PRODOM',
		       'smart'       => 'SMART',
		       'pfam'        => 'PFAM',
		       'ncbifam'     => 'NCBIFAMs',
		       'ncoils'      => 'COIL',
		       'seg'         => 'SEG',
		       'tmhmm'       => 'TMHMM',
		       'signalp'     => 'SIGNALP', 
		       'pirsf'       => 'PIRSF',
		       'superfamily' => 'SSF',
		       'gene3d'      => 'GENE3D',
		       'hmmpanther'  => 'PANTHER',
		       'hamap'       => 'HAMAP',
	       );

=head2 new

     Title    : new
     usage    : $m = new Domain2Interpro
     Function : constructor to provide the mapping data
     returns  : a Domain2Interpro object
     Args     : none

=cut

# constructor
sub new {
	my $class = shift;
	my %self;
        my ($ip_ids,$ip2name) = &get_ip_mappings();# hash of Databases hash of IDs
	$self{ip_ids}=$ip_ids;
	$self{ip2name}=$ip2name;
	bless \%self , $class;
}

=head2 get_method2database

     Title    : get_method2database
     usage    : $m->get_method2database('Tigrfam')
     Function : accessor to the class variable
     returns  : a string containind the database id
     Args     : method name

=cut

# mapping of method to database
sub get_method2database {
	my ($self,$key)=@_;
	return $method_database{lc($key)}}

=head2 get_mapping

     Title    : get_mapping
     usage    : $m->get_mapping(\%features)
     Function : default method to get mappings
     returns  : list containing arrayrefs of @hit = ( $ip_id, $start, $end, $hstart, $hend, $score, $evalue )
     Args     : lighweight feature: hashref->{method}->arrayreff(array($hid,$start,$end,$hstart,$hend,$score,$evalue))

=cut

# lets assume an hashref->{method}->arrayreff(array($hid,$start,$end,$hstart,$hend,$score,$evalue))
sub get_mapping{
	my ($self,$data)=@_;
	my @m=$self->get_motifs($data);
	return $self->merge_hits(@m);
}

=head2 get_motifs

     Title    : get_motifs
     usage    : $m->get_motifs
     Function : method to get mappings of motifs (unmerged)
     returns  : list containing arrayrefs of @hit = ( $ip_id, $start, $end, $hstart, $hend, $score, $evalue )
     Args     : lighweight feature: hashref->{method}->arrayreff(array($hid,$start,$end,$hstart,$hend,$score,$evalue))

=cut

# lets assume an hashref->{method}->arrayref(array($hid,$start,$end,$hstart,$hend,$score,$evalue))
sub get_motifs {
	my ($self,$motifs) =@_;
        # get the motifs
        my @_motifs;
        while( my ($method,$arefs)=each %$motifs) {
         foreach my $aref (@$arefs) {
           my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = @$aref;
           if ((lc($method) eq "pfam" ) && ( $hid =~ /(\w+)\.\d+/ )) { $hid = $1}
	   elsif (lc($method) eq "superfamily") { $hid = "SSF$hid"}
	   elsif (lc($method) eq "gene3d") {$hid = "G3DSA:$hid"}
           # convert Database ID to InterPro ID (if it is in InterPro)
           my $database = $method_database{lc($method)};
           if (exists $self->{ip_ids}{$database}{$hid}) {
             my $ip_id = $self->{ip_ids}{$database}{$hid};
             my @hit = ( $ip_id, $start, $end, $hstart, $hend, $score, $evalue );
             push @_motifs, [ @hit ];
           } 
          }
        }
	return @_motifs;
}

=head2 merge_hits

     Title    : merge_hits
     usage    : $m->merge_hits(@ip-features)
     Function : default method to get merge ip features
     returns  : list containing arrayrefs of @hit = ( $ip_id, $start, $end, $hstart, $hend, $score, $evalue )
     Args     : list containing arrayrefs of @hit = ( $ip_id, $start, $end, $hstart, $hend, $score, $evalue )

=cut

# merge overlapping hits with the same IP id
sub merge_hits {
	my ($self,@motifs)=@_;        
	my @merged;
        # sort the hits for this protein by ID and start position
        @motifs = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @motifs;
 
        # go through the hits for this protein looking for the same InterPro Id hits that overlap
        # then merge them
        my $prev_id = "";
        my ($merged_start,$merged_end, $merged_hstart, $merged_hend, $merged_score, $merged_evalue);
        my @merged_hit;
        foreach my $hit (@motifs) {
        my ( $ip_id, $start, $end, $hstart, $hend, $score, $evalue ) = @$hit;
        # is this a new ID or the same ID at a different part of the protein?
        if ($prev_id ne $ip_id || $start > $merged_end) {

          # save the merged hit
          if ($prev_id ne "") {
 	    @merged_hit = ( $prev_id, $merged_start, $merged_end, $merged_hstart, $merged_hend, $merged_score, $merged_evalue,$self->{ip2name}->{$prev_id});
	    push @merged, [ @merged_hit ];
          }

          # reset the values for the merged hit to be the values of the new ID
          $prev_id = $ip_id;
          $merged_start = $start;
          $merged_end = $end;
          $merged_hstart = $hstart;
          $merged_hend = $hend;
          $merged_score = $score;
          $merged_evalue = $evalue; 
        } else {
          # merge this hit into the merged hits
          # some results (e.g. prints) have zero hstart and hend, so overwrite these
          if ($merged_start > $start) {$merged_start = $start}
          if ($merged_end < $end) {$merged_end = $end}
          if ($hstart != 0 && ($merged_hstart == 0 || $merged_hstart > $hstart)) {$merged_hstart = $hstart}
          if ($merged_hend == 0 || $merged_hend < $hend) {$merged_hend = $hend}
          if ($merged_score < $score) {$merged_score = $score}
          if ($merged_evalue > $evalue) {$merged_evalue = $evalue}
        }
       }
       # save the last merged ID
       @merged_hit = ( $prev_id, $merged_start, $merged_end, $merged_hstart, $merged_hend, $merged_score, $merged_evalue ,$self->{ip2name}->{$prev_id} );
       push @merged, [ @merged_hit ];
       return @merged;
}

=head2 get_interpro

     Title    : get_interpro
     usage    : get_interpro('/tmp/')
     Function : get the interpro XML file and save it
     returns  : 
     Args     : directory to save the file in

=cut

#########################
# get the interpro file #
#########################

sub get_interpro {
   my $file = shift;
  `wget -O $file ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz`;
 }

=head2 get_ip_mappings

     Title    : get_ip_mappings
     usage    : get_ip_mappings
     Function : build the motif->ip mappings
     returns  : hash of database IDs , hash hashes of interpro to motif ids for each method
     Args     :

=cut


#########################################################
# reads in data for database ID to InterPro ID mapping  #
#########################################################
sub get_ip_mappings {
  # the interpro.xml file can be obtained from:
  # ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz

  my $file = "/tmp/interpro.$$.xml.gz";

  # get the interpro file from the EBI
  unlink "$file" if -e "$file";
  get_interpro($file);
 
  open (XML, "zcat $file|") || die "Failed to open pipe to zcat $file\n";
 
  my $in_member_list = 0;       # flag for in data ID section of XML file
  my ($IPid, $IPname, $this_db, $this_dbkey);
  my %ip_ids;       # hash of Databases hash of IDs
  my %ip2bla;       # ip_id2name

  while (my $line = <XML>) {
    if ($line =~ /<interpro id=\"(\S+)\"/) {
      $IPid = $1;
    } elsif ($line =~ m|<name>(.+)</name>|) {
      $ip2bla{$IPid}=$1;
    } elsif ($line =~ /<member_list>/) { # start of database ID section
      $in_member_list = 1;
    } elsif ($line =~ m|</member_list>|) { # end of database ID section
      $in_member_list = 0;
    } elsif ($in_member_list && $line =~ /<db_xref/) {
      ($this_db, $this_dbkey) = ($line =~ /db=\"(\S+)\" dbkey=\"(\S+)\" /);
      $ip_ids{$this_db}{$this_dbkey} = $IPid;
    }
  }
  close (XML);
  unlink $file;
  return (\%ip_ids,\%ip2bla);
}

=head2 DESTROY

     Title    : DESTROY
     usage    : $m->DESTROY or automatically called
     Function : cleanes up any tmp files
     returns  :
     Args     :

=cut


sub DESTROY {
  my $file ='/tmp/interpro.$$.xml.gz';
  unlink  $file if -e $file;
}

=head2 domain_lookup

     Title    : domain_lookup
     usage    : $m->domain_lookup('PFAM', $pf)
     Function : convert database entry to InterPro domain
     returns  : string or undef
     Args     : source db as string, db id as string.

=cut

sub domain_lookup {
    my $self = shift;
    my $db = shift;
    my $lookup = shift;
    return unless $lookup;

    my $domain = $self->{'ip_ids'}->{$db}->{$lookup};
    return $domain ? $domain : undef;
}

sub get_domain_name {
    my $self= shift;
    my $domain = shift;
    return $self->{'ip2name'}->{$domain} ? $self->{'ip2name'}->{$domain} : undef ;
}

1;
