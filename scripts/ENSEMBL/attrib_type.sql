-- MySQL dump 9.11
--
-- Host: ecs2    Database: caenorhabditis_elegans_core_39_150a
-- ------------------------------------------------------
-- Server version	4.1.12

--
-- Dumping data for table `attrib_type`
--

INSERT INTO `attrib_type` VALUES (1,'embl_acc','EMBL accession',NULL);
INSERT INTO `attrib_type` VALUES (2,'status','Status',NULL);
INSERT INTO `attrib_type` VALUES (3,'synonym','Synonym',NULL);
INSERT INTO `attrib_type` VALUES (4,'name','Name','Alternative/long name');
INSERT INTO `attrib_type` VALUES (5,'type','Type of feature',NULL);
INSERT INTO `attrib_type` VALUES (6,'toplevel','Top Level','Top Level Non-Redundant Sequence Region');
INSERT INTO `attrib_type` VALUES (7,'GeneCount','Gene Count','Total Number of Genes');
INSERT INTO `attrib_type` VALUES (8,'KnownGeneCount','Known Gene Count','Total Number of Known Genes');
INSERT INTO `attrib_type` VALUES (9,'PseudoGeneCount','PseudoGene Count','Total Number of PseudoGenes');
INSERT INTO `attrib_type` VALUES (10,'SNPCount','SNP Count','Total Number of SNPs');
INSERT INTO `attrib_type` VALUES (11,'codon_table','Codon Table','Alternate codon table');
INSERT INTO `attrib_type` VALUES (12,'_selenocysteine','Selenocysteine',NULL);
INSERT INTO `attrib_type` VALUES (13,'bacend','bacend',NULL);
INSERT INTO `attrib_type` VALUES (14,'htg','htg','High Throughput phase attribute');
INSERT INTO `attrib_type` VALUES (15,'miRNA','Micro RNA','Coordinates of the mature miRNA');
INSERT INTO `attrib_type` VALUES (16,'non_ref','Non Reference','Non Reference Sequence Region');
INSERT INTO `attrib_type` VALUES (17,'sanger_project','Sanger Project name',NULL);
INSERT INTO `attrib_type` VALUES (18,'clone_name','Clone name',NULL);
INSERT INTO `attrib_type` VALUES (19,'fish','FISH location',NULL);
INSERT INTO `attrib_type` VALUES (21,'org','Sequencing centre',NULL);
INSERT INTO `attrib_type` VALUES (22,'method','Method',NULL);
INSERT INTO `attrib_type` VALUES (23,'superctg','Super contig id',NULL);
INSERT INTO `attrib_type` VALUES (24,'inner_start','Max start value',NULL);
INSERT INTO `attrib_type` VALUES (25,'inner_end','Min end value',NULL);
INSERT INTO `attrib_type` VALUES (26,'state','Current state of clone',NULL);
INSERT INTO `attrib_type` VALUES (27,'organisation','Organisation sequencing clone',NULL);
INSERT INTO `attrib_type` VALUES (28,'seq_len','Accession length',NULL);
INSERT INTO `attrib_type` VALUES (29,'fp_size','FP size',NULL);
INSERT INTO `attrib_type` VALUES (30,'BACend_flag','BAC end flags',NULL);
INSERT INTO `attrib_type` VALUES (31,'fpc_clone_id','fpc clone',NULL);
INSERT INTO `attrib_type` VALUES (32,'KnownPCCount','protein_coding_KNOWN','Number of Known Protein Coding');
INSERT INTO `attrib_type` VALUES (33,'NovelPCCount','protein_coding_NOVEL','Number of Novel Protein Coding');
INSERT INTO `attrib_type` VALUES (34,'NovelPTCount','processed_transcript_NOVEL','Number of Novel Processed Transcripts');
INSERT INTO `attrib_type` VALUES (35,'PutPTCount','processed_transcript_PUTATIVE','Number of Putative Processed Transcripts');
INSERT INTO `attrib_type` VALUES (36,'PredPCCount','protein_coding_PREDICTED','Number of Predicted Protein Coding');
INSERT INTO `attrib_type` VALUES (37,'IgSegCount','total_Ig_segment_','Number of Ig Segments');
INSERT INTO `attrib_type` VALUES (38,'IgPsSegCount','Ig_pseudogene_segment_','Number of Ig Pseudogene Segments');
INSERT INTO `attrib_type` VALUES (39,'TotPsCount','total_pseudogene_','Number of Pseudogenes');
INSERT INTO `attrib_type` VALUES (40,'ProcPsCount','processed_pseudogene_','Number of Processed pseudogenes');
INSERT INTO `attrib_type` VALUES (41,'UnprocPsCount','unprocessed_pseudogene_','Number of Unprocessed pseudogenes');
INSERT INTO `attrib_type` VALUES (42,'KnwnPCProgCount','protein_coding_in_progress_KNOWN','Number of Known Protein coding genes in progress');
INSERT INTO `attrib_type` VALUES (43,'NovPCProgCount','protein_coding_in_progress_NOVEL','Number of Novel Protein coding genes in progress');
INSERT INTO `attrib_type` VALUES (44,'AnnotSeqLength','Annotated sequence length','Total length of annotated clone sequence');
INSERT INTO `attrib_type` VALUES (45,'TotCloneNum','Total number of clones','Total number of clones');
INSERT INTO `attrib_type` VALUES (46,'NumAnnotClone','Fully annotated clones','Number of fully annotated clones');
INSERT INTO `attrib_type` VALUES (47,'ack','Acknowledgement','');
INSERT INTO `attrib_type` VALUES (48,'htg_phase','High throughput phase','high throughput genomic sequencing phase');
INSERT INTO `attrib_type` VALUES (49,'description','Description','A general descriptive text attribute');
INSERT INTO `attrib_type` VALUES (50,'chromosome','Chromosome','chromosomal location for supercontigs that are not assembled');
INSERT INTO `attrib_type` VALUES (51,'nonsense','Nonsense Mutation','Strain specific nonesense mutation');
INSERT INTO `attrib_type` VALUES (52,'author','Author','Group resonsible for Vega annotation');
INSERT INTO `attrib_type` VALUES (53,'author_email','Author email address','Author email address');
INSERT INTO `attrib_type` VALUES (54,'remark','Remark','Annotation remark');
INSERT INTO `attrib_type` VALUES (55,'transcr_class','Transcript class','Transcript class');
INSERT INTO `attrib_type` VALUES (56,'KnownPTCount','processed_transcript_KNOWN','Number of Known Processed Transcripts');
INSERT INTO `attrib_type` VALUES (57,'ccds','CCDS','CCDS identifier');
INSERT INTO `attrib_type` VALUES (58,'initial_met','Initial methionine','Set first amino acid to methionine');
INSERT INTO `attrib_type` VALUES (59,'Frameshift  Fra','Frameshift modelled as intron',NULL);
INSERT INTO `attrib_type` VALUES (60,'NovelCDSCount','protein_coding_NOVEL','Total Number of Novel CDSs');
INSERT INTO `attrib_type` VALUES (61,'NovelTransCount','processed_transcript_NOVEL','Total Number of Novel transcripts');
INSERT INTO `attrib_type` VALUES (62,'PutTransCount','processed_transcript_PUTATIVE','Total Number of Putative transcripts');
INSERT INTO `attrib_type` VALUES (63,'PredTransCount','protein_coding_PREDICTED','Total Number of Predicted transcripts');
INSERT INTO `attrib_type` VALUES (64,'UnclassPsCount','pseudogene_NOVEL','Number of Unclassified pseudogenes');
INSERT INTO `attrib_type` VALUES (65,'KnwnprogCount','protein_coding_in_progress_KNOWN','Number of Known Genes in progress');
INSERT INTO `attrib_type` VALUES (66,'NovCDSprogCount','protein_coding_in_progress_NOVEL','Number of novel CDS in progress');
INSERT INTO `attrib_type` VALUES (67,'GeneNo_rRNA','rRNA Gene Count','Number of rRNA Genes');
INSERT INTO `attrib_type` VALUES (68,'GeneNo_knwCod','known protein_coding Gene Count','Number of known protein_coding Genes');
INSERT INTO `attrib_type` VALUES (69,'GeneNo_pseudo','pseudogene Gene Count','Number of pseudogene Genes');
INSERT INTO `attrib_type` VALUES (70,'Frameshift','Frameshift','Frameshift modelled as intron');

