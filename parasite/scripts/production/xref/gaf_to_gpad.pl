#!/usr/bin/env perl
# Switch the Uniprot GAF file into a GPAD file 
# Consumer: $ENSEMBL_CVS_ROOT_DIR/ensembl-production/modules/Bio/EnsEMBL/Production/Pipeline/GPAD/LoadFile.pm
# From ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README
#GAF2.1 files have the suffix .gaf and contain the following columns:
#
#	Column  Contents
#	1       DB
#	2       DB_Object_ID
#	3       DB_Object_Symbol
#	4       Qualifier
#	5       GO_ID
#	6       DB:Reference
#	7       Evidence Code
#	8       With (or) From
#	9       Aspect
#	10      DB_Object_Name
#	11      DB_Object_Synonym
#	12      DB_Object_Type
#	13      Taxon and Interacting taxon
#	14      Date
#	15      Assigned_By
#	16      Annotation_Extension
#	17      Gene_Product_Form_ID
#
#GPAD1.1 files have the suffix .gpa and contain the following columns:
#
#	Column  Contents
#	1       DB
#	2       DB_Object_ID
#	3       Qualifier
#	4       GO_ID
#	5       DB:Reference
#	6       Evidence Code
#	7       With (or) From
#	8       Interacting_Taxon_ID
#	9       Date
#	10      Assigned_By
#	11      Annotation_Extension
#	12      Annotation_Properties
# GAF 7 is not the same as GPAD 6, so we can't provide it, but we don't.
# Unused by Ensembl (in version 92):
#    GPAD 3, 5, 6, 7, 9, 11
# Ensembl gets nice files, with Ensembl transcript / protein ids assigned.
# We don't have that in the GAF files, but since this pipeline will be ran after xrefs, it's okay.

sub gaf_from_gpad {
  return join ("\t", 
   @_[1], 
   @_[2], 
   @_[4], 
   @_[5], 
   @_[6],
   "" , 
   @_[8],
   "" ,
   @_[14],
   @_[15],
   @_[16],
   sprintf("tgt_taxon=%s|go_evidence=%s\n", (@_[13] =~ s/taxon://r),@_[7] )
  );
}

while (<>){
  print join "\t", &gaf_from_gpad("", split( /\t/, $_, -1));
}
