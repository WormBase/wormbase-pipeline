#!/usr/local/bin/perl -w- -*-Perl-*-
#smlib2.pl of frequently used functions 
#validgene = no return unless a recognised gene or allele or rearrangement. returns the object
#validauthor checks that they have Surname INITS format or warn
#genetype checks if normal locus, rearr and returns the type of object
#validcalc for 2 point-crosses
#convertdate
#tidy_names removes fullsopts form naems, converts and/& into ,
#convertmapper
#reverseauthor takes in inits surname and reverses them.
#convert_author take sin potential authro and converts into ace format, uses other progs
#object_author
#convertmapper
#SINCE VALIDGENEARR AND VALIDGENEONE NEVER USED, HAVE TAKEN OUT

#######################################################################################
#VALIDATE C.ELEGANS GENE

#to check the format of nematode genes, separate file for now
#need to pass it either a single var or an array or  of mixed or pure data
#I can pass each cell of an array and report back on individual cells
#OR pass an array and report on whole array with any error messages.

#################################################################################################

sub validgenetwo #pared down, no report statements now, use validgeneone for these
{          
           local($g);
           $g = $_[0];
	   #print "this is doll-g $g\n\n";

	   if ($g =~ /^\w+\n/){ chop $g;  }    #remove return if there
           elsif ($g =~ /^\w+\s+/){ s/\s+//; }
          
           s/\(//g;                              #remove brackets from alleles
	   s/\)//g;
	   
  	if ($g =~/^[a-z][a-z][a-z]-\d+/)    #if C.elegans normal gene
	{	
                $g;	
	}
	elsif ($g  =~/^[a-z]+[DTC][pf]?\d+/) #if a rearrangemnt
        {
               $g;
	}
        elsif ($g =~/^[a-z]{1,2}P\d+/)        #if a polymorphism  DOES NOT WORK
	{	
               $g;
	}
        elsif ($g =~/^[a-z]{1,2}\d+/)        #if an allele
	{	
               $g;
	}
	#else no return unless a recognised gene or allele or rearrangement
}

######################################################################################

sub genetype
{
           #print "\nI'm in genetype function\n";
           local($g);
           $g = $_[0];

	   #print " (genetype) before chop = $g ";
	   if ($g =~ /\w+\n/){ chop $g;} #print "(genetype) after chop = $g";  }
	   
	   s/\(//g;
	   s/\)//g;
	   
  	if ($g =~/^[a-z][a-z][a-z]-\d+/)
	{	
                #print "Locus ";
		$classname = Locus;
	        $classname;          #return classname
	}
	elsif ($g  =~/^[a-z]+[DTC][pf]?\d+/) 
        {
               #print "Rearrangement ";
	       $classname = Rearrangement;
	       $classname;
	}
	#elsif ($g =~/^ayP1/)
        elsif ($g =~/^[a-z]{1,2}P\d+/)
	{	
               #print "Polymorphism ";
	       $classname = Polymorphism;
	       $classname;
	}
        elsif ($g =~/^[a-z]{1,2}\d+/)
	{	
               #print "Allele ";
	       $classname = Allele;
	       $classname;
	}
	#else 
        #{ 
        #     print  "gene $g  probably not C.elegans locus type  in line $.\n\n";	
        #}
}

####################################################################################

sub validcalc
{
  local($c);
  $c = "";
  $c = $_[0];
		#what on earth is wrong with this I put & instread of sub!

  if (( $c =~ /(Full|Backcross|Sex_all)/ ) && ($c =~ /\s*\d+\s+\d+\s+\d+\s+\d+\s*/))
  { 	print "calc is 4-part\n";
  	$c;
  }
	
  elsif ( ( $c =~ /(Recs_all)/ ) && ( $c =~ /\s*\d+\s+\d+\s+\d+\s*/ ) )
  {
  	print "calc is 3-part\n";
	$c;
  }

  elsif(($c =~/(One_all|One_rec|One_recombinant|One_let|Dom_let|Dom_semi|Dom_selected|Dom_one)/) 
		&& ( $c =~ /\s*\d+\s+\d+\s*/ ) )
  { 	print "calc is 2-part\n";
	$c;
  }

 elsif(($c =~ /Selected|Tested|Direct|Selected_trans|Back_one|Complex_mixed|Sex_one|Sex_cis/)
	&& ( $c =~ /\s*\d+\s+\d+\s*/ ) )
  { 	print "calc is 2-part\n";
	$c;
  }
#reports true or nothing at all put warnings in calling program
  #else
  #{ 
  #	print "something wrong in calc line\n";
  #}
}
	
###############################################################################################


#doesnot make allowance for american dates
#does not really deal with just the year
##################################################################################
sub convertdate	#takes an array (date) from main and converts to ace format
{	
	local (@date, $date);
	@date = @_;
	$len = @date;
	if ($len > 2) #days given as well as month and year
	{ 
		shift(@date); 
		print  @date,"\n";
	}
	@date[0,1] = @date[1,0];	#year-month now
	if ($date[0] =~ /\d+\s+/)
	{ 
		$date[0] =~ s/\s+//;
	}

	#if 2 digits for month carry on
	if ($date[1] =~ /[0-9]{2}/)
	{ 	
		print "Date $date[0]-$date[1]\n"; 
                $date = "$date[0]-$date[1]";
	}
	#elsif 1 digit, need to pad the month with 0
	elsif (@date =~ /\d+/)
	{	
		$date[1] = "0".$date[1];
		print "Date $date[0]-$date[1]\n";
                $date = "$date[0]-$date[1]";
	}
	else 
	{ 
		print "date wrong on line $.\n";
	}
	
	#else nothing to return
        $date;
}


	
###############################################################################################


#sub validauthor
#to take in author (object format) name(s)
# and check that they have Surname INITS format or warn
#do I need to be prepared to reverse name in this one
#NO :  each function should just do one job. (see sub object_author)
#could do french de ... or de la variation
#needs a bit of adjustment Bloch_Gallego WD and van der Biggs WH still a problem
###########################################################################################
	
sub validauthor
    {
		local($author) = "";
		
                print "validauthor doll0 = $_[0]\n";
		#this should not return Bloch-Gallego (no ints) but does so unless $anchor
		if ( ($_[0] =~ /^[A-Z]{1}[a-z]{1,}-[A-Z]{1}[a-z]{1,}\s+[A-Z]{1,4}\s*$/)
                || ($_[0] =~ /^[A-Z]{1}[a-z]{1,}-[a-z]{1,}\s+[A-Z]{1,4}\s*$/) )
		{		
			print "hyphenated surname $_[0]\n\n";
			$author = $_[0];
			$author;
		}

		elsif ($_[0] =~ /^(Von|von|van){1}\s*(der|den){0,1}\s*[A-Z]{1}[a-z]{1,}\s+[A-Z]{1,4}\s*$/)
		{   
			print "Dutch or German name $_[0] \n\n";
			$author = $_[0];
			$author;
		}

	        elsif ($_[0] =~ /^(Mac|Mc)[A-Z]{1}[a-z]{1,}\s+[A-Z]{1,4}\s*$/)
	        {
			print "Scottish author $_[0]\n\n";
			$author = $_[0];
			$author; 
		}

 		elsif ($_[0] =~ /^[A-Z]{1}[a-z]{1,}\s+[A-Z]{1,4}\s*$/)
		{
			print "normal author $_[0]\n\n";
			$author = $_[0];
			$author; #this return occurs only  if I've added $author = $_[0]
		}
                elsif ($_[0] =~ /\s*[A-Z]{1}[a-z]{1,}\s+[A-Z]{1,4}\s*/)
	      #authors, which have been split, don't respond to ^ and $
		{
			print "probably normal author $_[0]\n\n";
			$author = $_[0];
			$author; #this return occurs only  if I've added $author = $_[0]
		}
		#else
		#{
		#	print "irregular name $_[0] on line $.\n\n"; #or could warn...
		#}
}#WARN IN  MAIN PROGRAM
#END SUB VALIDAUTHOR
#################################################################################	

sub tidy_names
{
  
        local ($fullname) = ""; 
        $fullname = $_[0];
        #print "fullname IN tidy_author SUB = $fullname\n";
	$fullname =~ s/\.//g;# && print "\nIN tidy_author: removed fullstops from $fullname\n";

       		#if more than 1 author
        $fullname =~ /\band\b/ && $fullname =~ s/and/,/;# && print "changed and in $fullname\n"; 
        $fullname =~ /&/ && $fullname =~ s/&/,/; # && print "changed & in $fullname\n";
        $fullname =~ /\// && $fullname =~ s/\//,/; #&& print "changed / in $fullname\n";

	$fullname;
}#END SUB TIDY_NAMES

################################################################################

#to take in name(s) and convert them into author object names
#using other sub-routines when necessary

sub convert_author
{
 local($author,$fullname,$part,$objname,$firstname,$len,$objname,$objname2,$firstname,$init1,$name_type) = "";
 local (@name) = "";   

   $fullname = $_[0];
   #print "\n\nfullname IN convert_author SUB = $fullname\n";
   $name_type = $_[1];
   #print "name_type IN convert_author SUB = $name_type\n";

#remove unnecessary punctuation, convert connectors between authors to ,
        if($fullname =~ /[\.|and|&|\/|]/) 
        {                      # &tidy_names RETURNS  a name
              $fullname = &tidy_names($fullname)
              || warn " something gone wrong with $fullname\n";
        }

#split into individual authors.
# if there are more than 1 ##NEED TO ATTACH SURNAME TO BOTH NAMES ..NOT DONE
        if($fullname =~ /,/)
        {
               print "double author = $fullname\n";
               @name = split(/,/,$fullname);
	       foreach $part (0..$#name)
	       {
			chop if ($name[$part] =~/\n/);
                        $name[$part] =~ s/^\s+//g;	
                        #$name[$part] .= "\n"; this goes wrong later
			print "name[$part]= $name[$part]\n";
	       }
	       $len = @name;
	       $len >2 && print "\nnumber of names = $len\n"; 
	       $len >3 && print "trouble possible with $fullname\n"; 
	}
        else
        {
             $name[0] = $fullname;
        }
	
        foreach $n (0..$#name)
        {
        #If in shortened name form (INTS Surname) does it need reversing? 
                     #&reverseauthor RETURNS a name
	    $name[$n] = &reverseauthor($name[$n]);
            #print "From SUB reversed_author now = $name[$n]\n";

                            #&object_author RETURNS a name
	    #does it need converting from full to abbreviated name
	    #if consists of Name INITS, doesn't need objectifying
            #if($name[$n] =~ /^[A-Z]{1}[a-z]+\s+[A-Z]{1,3}\s*/)
#at least two names = like firstname surname or firstname middlename  

            if(($name[$n] =~ /^[A-Z]{1}[a-z]+\s+[A-Z]{1}[a-z]+/)
            || ($name[$n] =~ /^[A-Z]{1}[a-z]+\s+[A-Z]{1}\s+\w+/)  #does not allow for van/der/den/de as 2nd name
            || ($name[$n] =~ /^[A-Z]{1}[a-z]+-[a-z]+\s+[A-Z]{1}\w+/) # hyphenated first name
            || ($name[$n] =~ /^[A-Z]{1}[a-z]+-[A-Z]{1}[a-z]+\s+[A-Z]{1}\w+/) # hyphenated first name with 2 CAPS
	    || ($name[$n] =~ /^[A-Z]{1}[a-z]+\s+[de|van|der|den]{1,}\s+\S+/) 
            || ($name[$n] =~ /^[A-Z]{1}\s+[A-Z]{1}[a-z]+\s+\w+/) # R John Lye
            || ($name[$n] =~ /^[A-Z]{1}[a-z]+\s+O\'\w+/) ) # Bernard O'Connell
#or name = initial surname
            {
	         #print "probably needs full-->object name doing\n";
                 $name[$n] = &object_author($name[$n]);
                 print "objectified author = $name[$n]\n";
             }
	    else
            {
	        print "CARE: standard or WIERDO $name[$n] IN SUB CONVERT AUTHOR\n";
            }
            
            #validate the name
                            #&validauthor RETURNS  a name
	    $validname = &validauthor($name[$n]);
            print "validatedauthor = $validname\n"; 
            if ($validname =~ /\n/) {chop $validname; }
            #needs to send back a print statement not a var. 
            #more than one may need to go back
	    print "$name_type ";
            print "\"$validname\"\n"; 
            #this gives \n to single names on top of what is there
        }
}#END SUB CONVERT AUTHOR

####################################################################

#implies that name is already abbreviated to....
#INPUT = Author : "Initials Surname"
#some are double author but should have been resolved earlier

#########################################################################
sub  reverseauthor
{
   local($inits) = ""; local($surname) = "";
   local($fullname) = ""; local($author) = "";
   local($middlename) = "";
   $fullname = $_[0];
   #print "\nIN SUB reverseauthor fullname = $fullname\n";
   $fullname =~ s/\"//g; #remove any quotes
   $fullname =~ s/^([A-Z]{1})\s+([A-Z]{1})\s+([A-Z]{1})(\s+\.+)/$1$2$3$4/;
   $fullname =~ s/^([A-Z]{1})\s+([A-Z]{1})(\s+\w+.+)/$1$2$3/;
   if($fullname =~ /^([A-Z]{1,3})\s+(\w+)/)
   {
       print "IN SUB reverseauthor: author with just inits & surname to reverse\n";
       
       $inits = $1; $surname = $2;
   }
   elsif($fullname =~ /^([A-Z]{1}\w+)\s+(\w+)\s+(\w+)/ )
   {
          print "IN SUB reverseauthor : AUTHOR firstname and 2 surnames\n";
          $inits = $1; $middlename = $2; $surname = $3; 
   }
   elsif($fullname =~ /^([A-Z]{1}\w+)\s+(\w+)$/ )
   {
          print "IN SUB reverseauthor : AUTHOR firstname and surname\n";
   }
   
   else
   {
          print "IN SUB reverseauthor : AUTHOR VERY FUNNY\n";
   }

   if (($inits =~ /\w+/) && ($surname =~ /\w+/) && ($middlename !~ /\w+/) )
   {

        print "initials are $inits\t";
        print "surname is $surname\n";
        $author = $surname." ".$inits;
        print "\nAUTHOR NOW =",$author,"\n";
        #put through author validation with full warnings
         &validauthor($author) && print "AUTHOR= $author\n\n\n"
	 || warn "funny author line $. \n ";

         $author;
  }
  else
  {
        print "returning original full_name\n";
        $fullname;
  }
}#END SUB REVERSEAUTHOR
######################################################################

###########SUB-ROUTINE FOR REVERSING FULL AUTHOR NAMES###################
###########AND CONVERTING THEM INTO OBJECT NAMES#########################

#sub object_author mainly used for subscribers to WBG
#make object name by reversing surname and firstnames and 
#chopping off first name to its initial
#then needs running through validateauthor routine

#fullnames formats recognised by progam : Fred Bloggs, Fred M Bloggs, Fred  La Bloggs, 
#Fred M. Bloggs, Fred Martin Bloggs, Fred von Bloggs,  Fred Martin-Bloggs, Fred de Bloggs;
#not recognised yet : Fred van der Bloggs; Fred de la Bloggs; i.e. 4 part names

#if firstname is CAPfullstop and second name is a fullname its probably a givenname.
#if firstname is Capletters, second name is Capletters, third is surname
#then second name assumed to be WHAT????
#Juan E. Abrahante Llorens NOT coped with except by a warning and giving no full_name



sub object_author
{
    local($author) = ""; local ($fullname) = ""; 
    local (@name) = ""; local ($part) = ""; local ($len) = "";
    local ($objname) = ""; local ($objname2) = ""; 
    local ($firstname) = ""; local ($init1) = "";
    local ($name1,$name2,$name3) = "";
  if($_[0]=~/^[A-Z]{1}/)#AT LEAST ITS NOT BLANK!
  {
	$fullname = $_[0]; 
	#print "full_name in SUB object_author = $fullname\n";#this is already done in main
	
#split name into parts, keep last name as surname = object first part	
	@name = split(/ /,$fullname);
	foreach $part (0..$#name)
	{
		print "name[$part]= $name[$part]\n";
	}
	$len = @name;
        if ($len ==2) { print "\nlength of name = TWO\n";}  
	if ($len >2) { print "\nlength of name = $len\n";} 
	if ($len >3) { print "serious trouble possible in $fullname, line $.\n"; }
	$objname = $name[$len-1]; #i.e. surname
	#print "objname initially = $objname\n";

#take first name in author as given name. take only first letter.
#should work for Sylvia /  S. variants
#@init = split (//,$firstname);$init2 = shift(@init);print "init2= $init2\t...
		
	if ($name[0] !~ /-/)
        {
           $firstname  = $name[0];
	   $init1 = substr($firstname,0,1);
        }
        elsif($name[0] =~ /-/) 
        {
            $firstname  = $name[0];
            $init1 = substr($firstname,0,1);
            $firstname =~ /\w+-(\w+)/;
            $extra = $1;
            $init1 .= "-";
            $init1 .=substr($extra,0,1);
        }
 
	#print  "init1= $init1\n";
	$objname .= " "; 
	$objname .= $init1; #the basic object name structure for all types
        #print "objname plus initial = $objname\n";
#then deal with rest of name if it exists, i.e extra initials, extra parts to surname.
        if($len ==2) {$objname;}

#if middlename = initial and has already had fullstop removed
	elsif(($len>2)&&($name[1]=~/([A-Z]){1}\b\s*/))
        { 
                  $name1 = $objname;
	          $name1 .= $name[1];
                  #print  "middle name = initial\n"; 
                  $name1;
	}

#if a full middle name, could be second given name, woman's previous name, part of complex surname. 
#assume its part of surname and treat as de von van etc.
                
	elsif(($len >2) && ($name[1] =~ /[\w]{2,}/))
	{	
	     print "non-standard name, $fullname, line $.\n";
	     print "a middle name $name[1] in full\n";
	     print "before= $objname\n";
			
#having trouble with ORs here, so copying syntax from validauthor
#/^(Von|von|van){1}\s*(der|den){0,1}\s*[A-Z]{1}[a-z]{1,}\s+[A-Z]{1,4}\s*$/)

		if($name[1] =~ /(van|Van|La|Las|De|de|von|Von|Dal){1}\s*(den|der|La){0,1}/)
		{
		      $name2 = $name[1]." ".$objname;
                      #print "final continental object name in SUB-ROUTINE = $name2\n\n";
	              $name2;
                }
                else #assuming its a given middle name
                {
                      $init2 = substr($name[1],0,1);
		      $name3 = $objname.$init2;
                      #print "final 'dunno' object name in SUB-ROUTINE object-author = $objname\n\n";
		      $name3;
                }
         }
#last print only really needed for funny names, hence included here.
#must modify to allow for van/von der/den, de la.
         elsif(($len > 3) && ($name[1] =~ /[\w]{2,}/))
         {
            #print "very funny name based on $objname\n";
	    $objname;
         }
		
    }#endif there is a name present
 
}#end SUB

#########################################################################################

#there are infinite variations in the mapper name presentation
# one objectified name Bloggs SB.... THE NORM
# Names separated by  : and or &
# Harris, JM
#J.Weems/D.Shakes
#2 names separated by comma is taken to be the norm and othe variations converted to this

#although sending back print messages is simple, it dependes on whether file handles
#are being used in the calling program so makes ithem interdependent.
#just returning an author is better.
#however there is a problem with two authors
#since a sub should only do one job, it should not split authors
#and check them

sub convertmapper	#takes a VAR which could be split later into array
{
        print "from convertmapper program\n";
        local $author;
	local @author;
	$author = $_[0];
        print "author/mapper = $author\n";
         #if($author !~ /[A-Z{1}[a-z]+\s+[A-Z]{1,}$/) { warn "could be author problem\n"; }
        $author =~ s/\.//g;

			#if more than 1 author
        $author =~ /and/ && $author =~ s/and/,/ && print "boo\t$author\n"; 
        $author =~ /&/ && $author =~ s/&/,/ && print "boo\t$author\n";
        $author =~ /\// && $author =~ s/\//,/ && print "boo\t$author\n";

	if ($author =~ /,/) #could just have , here
	{ 
		@author = split (/, /,$author);
		foreach $m (0..$#author)
		{	
	                print  $author[$m],"\n";
			&validauthor($author[$m]) && print OUTFILE "Mapper \"$author[$m]\"\n" 
			|| print "bad author $author[$m] in $history\n";
		}
	}
			#if only 1 author
	elsif ($author =~ /\w+[^,]/) 
	{ 
		if (&validauthor($author))
                {
		    $author;
                }
                else
                {
                    print "bad author $author\n";
                }
                  # && print OUTFILE "Mapper \"$author\"\n" 
		#|| print "bad author $author in $history\n";
	}
			#if no author
	elsif ($author =~ /[^w+]/)
	{
		print OUTFILE "-C \"mapper missing\"\n";
	}
			
}
#########################################################################################

#########################################################################################
###############################################################################################
###############################################################################
1;

