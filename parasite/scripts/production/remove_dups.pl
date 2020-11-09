use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Data::Dumper;

my (%descriptions, %new_descriptions, $host, $user, $db, $pass, $port, $updatedb);

GetOptions(
	'host=s'		=> \$host,
	'dbname=s'		=> \$db,
	'pass=s'		=> \$pass,
	'port=i'		=> \$port,
	'user=s'		=> \$user,
	'updatedb'		=> \$updatedb,
) || die ("check command line params\n");


my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $host,
    -user   => $user,
    -dbname => $db,
    -pass   => $pass,
    -port   => $port,
      );

#get gene descriptions
my $dump_sth = $dba->dbc->prepare("SELECT gene_id, description FROM gene WHERE description IS NOT NULL");
$dump_sth->execute;
while( my ($gene_id, $orig_description) = $dump_sth->fetchrow_array){
	$descriptions{$gene_id} = $orig_description;
}
$dump_sth->finish;

#remove duplicate sections
foreach my $gene_id (keys %descriptions){
	my $source;
        if ($descriptions{$gene_id} =~ /(\s*\[Source:.*\])$/){
                $source = $1;
        }
	else{ die "couldn't parse source from $descriptions{$gene_id}"; }

        (my $orig_desc = $descriptions{$gene_id}) =~ s/(\s*\[Source:.*\])$//;
        my @orig_descs = split (";", $orig_desc);

        #remove leading and trailing white space
        @orig_descs = map  { (my $this_desc = $_ ) =~ s/\s*$//; $this_desc } @orig_descs;
       	@orig_descs = map  { (my $this_desc = $_ ) =~ s/^\s*//; $this_desc } @orig_descs;
        
        my @unique_descs = uniq(@orig_descs);
        my $new_desc = join("; ", @unique_descs)."$source";
	if ($new_desc ne $descriptions{$gene_id}){
		$new_descriptions{$gene_id} = $new_desc;
	}	
}

#if we have new descriptions, log them in a file
#also write to database if flag is given
if (%new_descriptions){
	open(my $fh, '>', "$db.edited_descriptions.txt");
	my $update_sth = $dba->dbc->prepare("UPDATE gene SET description = ? WHERE gene_id = ?");
	foreach my $gene_id (keys %new_descriptions){
		print $fh "$gene_id\t$descriptions{$gene_id}\t$new_descriptions{$gene_id}\n";
		if ($updatedb){
			$update_sth->execute($new_descriptions{$gene_id}, $gene_id);
			$update_sth->finish;
			print "Updated $gene_id to $new_descriptions{$gene_id}\n";
		}
	}
}

