#!/software/bin/perl -w
use strict;
use lib $ENV{'CVS_DIR'};
#use NameDB_handler;
use Getopt::Long;
use Log_files;
use Ace;
use Wormbase;

=pod

=head batch_parse_allele_form.pl

=item Options:

  -file file containing a concatination of all the emails from the form <Mandatory>

    Simple FORMAT for reading in here, but just concatinate the emails as the script pulls in the data based on the html in the raw file:
    Your Name	John Smith -- WBPerson1234567
    Your E-mail Address	john@smith.uk
    PubMed ID	123456
    Allele Name	abc1
    Gene Name	test-1 -- WBGene12345678
    Sequence Name	clone.number
    Type of Alteration	Point Mutation
    Alteration Details	a to a
    Type of Mutation	Bad, real bad
    30 bp upstream	tttttttttttttttttttttttttttttt
    30 bp downstream	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

    Your Name	John Smith -- WBPerson1234567
    Your E-mail Address	john@smith.uk
    PubMed ID	1234567
    Allele Name	abc2
    Gene Name	test-2 -- WBGene12345679
    Sequence Name	clone.number
    Type of Alteration	Point Mutation
    Alteration Details	a to a
    Type of Mutation	Bad, real bad
    30 bp upstream	tttttttttttttttttttttttttttttt
    30 bp downstream	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
  
  -domain    set domain to whatever (Gene Variation Feature)
  -debug     limits to specified user <Optional>
  -species   can be used to specify non elegans genes
  -load      loads the resulting .ace file into geneace.
  -test      use the test nameserver
  -ns        Kill's the gene in the nameserver as well as producing the .ace 
             file for geneace
  -user      username                 <Manditory if using -ns>
  -password  password                 <Manditory if using -ns>

e.g. perl batch_parse_allele_form.pl -file deathrow.txt [simple example]

     perl batch_parse_allele_form.pl -u fred -p secret -file deathrow.txt -ns -test -debug mt3

=cut

my ($USER,$PASS, $test, $file, $ns, $debug, $load, $transposon,$seq_name,$wormbase,$out);
my $domain='Gene';
my $species = 'elegans';
GetOptions(
    'user:s'     => \$USER,
    'password:s' => \$PASS,
    'test'       => \$test,
    'file:s'     => \$file,
    'outfile:s'  => \$out,
    'ns'         => \$ns,
    'debug:s'    => \$debug,
    'load'       => \$load,
    'domain:s'   => \$domain,
    'species:s'  => \$species,
    ) or die;

my $log;
if (defined $USER) {$log = Log_files->make_log("NAMEDB:$file", $USER);}
elsif (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}
my $DB;
my $db;
my $ecount;

$wormbase = Wormbase->new("-organism" =>$species, -debug => $debug, -test => $test);
my $database = "/nfs/production/panda/ensemblgenomes/wormbase/DATABASES/geneace";
$log->write_to("TEST mode is ON!\n\n") if $test;
my $tace            = $wormbase->tace;        # TACE PATH

##############################
# warn/notify on use of -load.
##############################
if (!defined$load) {$log->write_to("You have decided not to automatically load the output of this script\n\n");}
elsif (defined$load) { $log->write_to("Output has been scheduled for auto-loading.\n\n");}

my $ace = Ace->connect (-path => $database,
                       -program => $tace) || die "cannot connect to database at $database\n";
my $output;
if ($out){
    $output = $out;
}
else {
    $output = $file.".out";
}

#open file and read
open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">$output") or $log->log_and_die("cant write output: $!\n");
my($person,$remark,$tflag,);
my $count;

my ($name, $wbperson, $public_name, $varname,$gene, $wbgene, $seq, $clone, $type_alt, $type_mut, $alt_det, $mut_det, $strain, $flank1, $flank2,$method,$geno,$mutagen,$forward,$comment,$pubmed,$obj_method,$raw_name);
while(<FILE>){
    chomp;
    
    if (/Your Name<\/td><td>(.+)<\/td>$/) { 
        $name = $1;
        #print "Name = $name\n"; 
        if ($name =~ /(WBPerson\d+)/) {
            $wbperson = $1;
        #    print "WBPerson ID = $wbperson\n";
            next;
        }
        else {
            next;
        }
    }
    elsif (/Your E-mail Address<\/td><td>(.+)<\/td>/) {
        #print "Email = $1\n";
        next;
    }

    elsif (/PubMed ID<\/td><td>(.+)<\/td>/) {
        #print "Pubmed_ID = $1\n";
        $pubmed = $1;
        next;
    }
    elsif (/Allele Name<\/td><td>(.+)<\/td>/) {
        #print "Allele name = $1\n";
        $raw_name = $1;
        
        if ($raw_name =~ /(\w+\d+)\s+\-\-\s+(WBVar\d+)/) {
            $public_name = $1;
            $varname = $2;
            next;
        }
        elsif ($raw_name =~ /(\w+\d+$)/) {
            $public_name = $1;
            next;
        }
        next;
    }
    elsif (/Gene Name<\/td><td>(\S+)<\/td>/) {
        #print "Gene name = $1\n";
        $gene = $1;
        if (/(WBGene{8})/) {
            $wbgene = $1;
        }
        next;
    }
    elsif (/Sequence Name<\/td><td>(\S+)<\/td>/) {
        #print "Seq name = $1\n";
        $seq = $1;
        if ($seq =~ /(\S+)\.\d+/) {
            $clone = $1;
        }
        next;
    }
    elsif (/Type of Alteration<\/td><td>(.+)<\/td>/) {
        #print "Type = $1\n";
        $type_alt = $1;
        next;
    }
    elsif (/Type of Mutation<\/td><td>(.+)<\/td>/) {
        #print "Type_Mut = $1\n";
        $type_mut = $1;
        next;
    }
    elsif (/Alteration Details<\/td><td>(.+)<\/td>/) {
        #print "Alt_type = $1\n";
        $alt_det = $1;
        next;
    }
    elsif (/Mutation Details<\/td><td>(.+)<\/td>/) {
        #print "Mut detail = $1\n";
        $mut_det = $1;
        next;
    }
    elsif (/30 bp upstream<\/td><td>(.+)<\/td>/) {
        #print "Flank_1 = $1\n";
        $flank1 = lc($1);
        next;
    }
    elsif (/30 bp downstream<\/td><td>(.+)<\/td>/) {
        #print "Flank_2 = $1\n";
        $flank2 = lc($1);
        next;
    }
    elsif (/Strain<\/td><td>(.+)<\/td>/) {
        #print "Strain = $1\n";
        my @strain_id_obj = $ace->fetch(-query => "FIND Strain where Public_name = $strain");
        if (defined $strain_id_obj[0]) {
            my $strain_id = $strain_id_obj[0]->name;
            $strain = $strain_id;
        }
        else {
            $strain = $1;
        }
        next;
    }
    elsif (/Production Method<\/td><td>(.+)<\/td>/) {
        #print "Prod_method = $1\n";
        $method = $1;
        next;
    }

    elsif (/Genotype<\/td><td>(.+)<\/td>/) {
        #print "Genotype = $1\n";
        $geno = $1;
        next;
    }
    elsif (/Mutagen<\/td><td>(.+)<\/td>$/) {
        #print "Mutagen = $1\n";
        $mutagen = $1;
        next;
    }
    elsif (/Forward genetics<\/td><td>(.+)<\/td>/) {
        #print "Forward_genetics = $1\n";
        $forward = $1;
        next;
    }
    elsif (/Comment<\/td><td>(.+)/) {
        #print "Comment = $1\n";
        $comment = $1;
        next;
    }
    elsif (/<\/td><td>(.+)<\/td>$/) { 
        #print "Caught this = $1\n";
        next;
    }
    elsif (/<\/table><br\/><br\/>/){    
        #print "Record_end\n\n"; 
    }
    else {
        #print "Ignoring $_\n";
        next;
    }
    

# my (, $type_alt, $type_mut, $alt_det, $mut_det
    if (/<\/table>/) {
        if (defined $varname){
            print ACE "Variation : \"$varname\"\n";
            my $status_obj = $ace->fetch(Variation=>$varname);
            my $seqstatus = $status_obj->SeqStatus->name;
            if ($seqstatus =~ /Sequenced/) {
                print "WARNING - $varname is already curated\n";
            }
        }
        else {
            
# Query geneace for WBVar name, if one exists return it, if one doesnt put in a holding ID (Not Implemented)

            my $obj = $ace->fetch(Variation_name=>$public_name);
            if (defined $obj) {
                my $wbvar_id = $obj->Public_name_for->name; 
                print ACE "Variation : \"$wbvar_id\"\n";
            }
            else {
                print ACE "Variation : \"$public_name\"\n";
            }
        }
        if (defined $wbperson) {print ACE "Evidence Person_evidence $wbperson\nPerson $wbperson\n";}
        print ACE "Public_name \"$public_name\"\nSpecies \"Caenorhabditis elegans\"\n";
        if ($clone) {
            print ACE "Mapping_target $clone\n";
        }
        else {
            print ACE "Mapping_target $seq\n";
        }
        if (defined $flank1){
            print ACE "Flanking_sequences \"$flank1\" \"$flank2\"\n";
            print ACE "SeqStatus Sequenced\n";
        }
        if (defined $wbgene) {
            if (defined $type_mut) {
                print ACE "Gene $wbgene $type_mut\n";
            }
            else {
                print ACE "Gene $wbgene\n";            
            }
        }

        if (defined $type_alt){
            if ($type_alt =~ /Engineered Allele/){
                #print "Ignoring type_alt Engineered Allele\'n" if ($debug);
            }
            elsif ($type_alt =~ /Point Mutation/){
                print ACE "Type_of_mutation Substitution $alt_det\n";
            }
            elsif ($type_alt =~ /Deletion/){
                print ACE "Type_of_mutation Deletion\n";
            }
            elsif ($type_alt =~ /Insertion/){
                print ACE "Type_of_mutation Insertion\n";
            }
        }
#            Insertion + Deletion
        if (defined $gene) {
            print ACE "Gene $gene\n" unless ($wbgene); 
        }
        if (defined $strain) {
            print ACE "Strain $strain\n";
        }
        if (defined $comment) {
            print ACE "Remark \"$comment\"\n";
        }
        if (defined $forward) {
            print ACE "Forward_genetics \"$forward\n";
        }
        if (defined $mutagen){
            print ACE "Mutagen \"$mutagen\"\n";
        }
        if (defined $method) {
            print ACE "Production_method $method\n";
        }

        if (defined $pubmed) {
#try and find the WBPaper ID
            my @paper_obj = $ace->fetch(-query=>"find Paper where Database AND NEXT AND NEXT AND NEXT = $pubmed");
            if (@paper_obj) {
            my $WBpaper = $paper_obj[0]->name;
            print ACE "Reference $WBpaper\n";
            }
            else {
                print "\/\/Reference $pubmed\n";
            }
        }
        #try and guess the method
        if ($method =~ /CRISPR_Cas9/) {
            $obj_method = "Engineered_allele";
        }
        elsif ($comment =~ /CRISPR/) {
            $obj_method = "Engineered_allele";
        }
        elsif ($type_alt =~ /Point Mutation/){
            $obj_method = "Substitution_allele";
        }
        elsif ($type_alt =~ /Substitution/) {
            $obj_method = "Substitution_allele";
        }
        elsif ($type_alt =~ /Deletion/){
            $obj_method = "Deletion_allele";
        }
        elsif ($type_alt =~ /Insertion/){
            $obj_method = "Insertion_allele";
        }
        elsif ($type_alt =~ /Substitution/){
            $obj_method = "Substitution_allele";
        }
        else {
            $obj_method = "Allele"; 
        }
        print ACE "Remark \"alt_det = $alt_det mut_det = $mut_det\"\n";
#standard tags
        print ACE "Sequenced\nSpecies \"Caenorhabditis elegans\"\nLive\n";
        print ACE "Method $obj_method";
        print ACE "\n\n";
        

        if (defined $geno) {
            if (defined $strain) {
                # Query geneace for a WBStrain name, if one exists return it, if one doesn't exist request a new one. (Not Implemented)
		# Work out the Lab from the Strain_name
		my $lab; 
		if ($strain =~ /([A-Z]+)\d+/){
		    $lab = $1;
		}
		else {undef $lab;}
                my @strain_obj = $ace->fetch(-query=>'find Strain where Public_name = $strain');
                if (@strain_obj){
                    my $strainID = $strain_obj[0]->name;
                    if (defined $strainID) {
                        print ACE "Strain : \"$strainID\"\nPublic_name $strain\nGenotype \"$geno\"\nSpecies \"Caenorhabditis elegans\"\n";
			if (defined $lab){print ACE "Location $lab\n\n";}
			else {print ACE "\n";}
                    }
                }
                else {
                    print ACE "Strain : \"$strain\"\nPublic_name $strain\nGenotype \"$geno\"\nSpecies \"Caenorhabditis elegans\"\n";
		    if (defined $lab){print ACE "Location $lab\n\n";}
		    else {print ACE "\n";}
                } 
            }
            else {
                print "\n\n\\\\No Strain for genotype exit without printing\n\n";
            }        
        }
        ($name, $wbperson, $public_name, $varname,$gene, $wbgene, $seq, $clone, $type_alt, $type_mut, $alt_det, $mut_det, $strain, $flank1, $flank2,$method,$geno,$mutagen,$forward,$comment,$pubmed,$obj_method,$raw_name) = ();
    }
}
        
&load_data if ($load);
$log->write_to("5) Check $output file and load into geneace.\n") unless ($load);
$log->mail();

sub load_data {
  # load information to $database if -load is specified
  $wormbase->load_to_database("$database", "$output", 'batch_parse_allele_form.pl', $log, undef, 1);
  $log->write_to("5) Loaded $output into $database\n\n");
  print "Finished!!!!\n";
}
