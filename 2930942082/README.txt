
Inside each <taxon_oid>.tar.gz bundle file:

<taxon_oid>.fna - FASTA nucleic acid file for taxon.
<taxon_oid>.genes.fna - FASTA nucleic acid file for genes.
<taxon_oid>.genes.faa - FASTA amnio acid file for genes.
<taxon_oid>.intergenic.fna - FASTA for intergenic regions.
<taxon_oid>.gff - Tab delimited in mostly GFF3 format for genes.
  (Strict conformance is not guarantted, esp. in type and attributes fields.)
<taxon_oid>.cog.tab.txt - Tab delimited file for COG annotation.
<taxon_oid>.kog.tab.txt - Tab delimited file for KOG annotation.
<taxon_oid>.pfam.tab.txt - Tab delimited file for Pfam annotation.
<taxon_oid>.tigrfam.tab.txt - Tab delimited file for TIGRFAM annotation.
<taxon_oid>.ipr.tab.txt - 
  Tab delimited file for "other" (Non-Pfam/TIGRFAM) InterPro hits.
  (Optional)
<taxon_oid>.ko.txt - Tab delimited file for KO and EC annotation.
<taxon_oid>.signalp.txt - Tab delimited file for signal peptide annotation.
<taxon_oid>.tmhmm.txt - Tab delimited file for transmembrane helices.
<taxon_oid>.xref.tab.txt - Tab delimited file for external references.
     (Data is spotty and optional.)
<taxon_oid>.crispr.txt - Tab delimited file for CRISPR details.

(For some of the smaller genomes, e.g., viruses, not all
 annotation files is present.  If the genome has no
 annotation of a certain type, e.g. no TIGRFAM's, 
 the annotation file <taxon_oid>.tigrmfam.tab.txt is not there.
 Some annotations are not done at all, e.g. InterPro for metagenomes.
 Any file that does not have annotations or has no data will
 not be present.)

------------
Structure of each tab delimited file:

<taxon_oid>.gff
-- seqid - Sequence ID
-- source - version of IMG database
-- type - feature type
-- start_coord - starting coordinate
-- end_coord - ending coordinate
-- score - NA
-- strand
-- phase - NA
-- attributes - ID=<gene_oid>;locus_tag=<locus_tag>;product=<product name>

<taxon_oid>.cog.tab.txt (from NCBI RPSBLAST)
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- percent_identity - Perceent identity of aligned amino acid residues
-- query_start - Start coordinate of alignment on query gene
-- query_end - End coordinate of alignment on query gene
-- subj_start - Start coordinate of alignment on subject sequence
-- subj_end - End coordinate of alignment on subject sequence
-- evalue - Expectation value
-- bit_score - Bit score of alignment
-- cog_id - COG identifier
-- cog_name - COG name
-- cog_length - Length of COG consensus sequence

<taxon_oid>.kog.tab.txt (from NCBI RPSBLAST)
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- percent_identity - Perceent identity of aligned amino acid residues
-- query_start - Start coordinate of alignment on query gene
-- query_end - End coordinate of alignment on query gene
-- subj_start - Start coordinate of alignment on subject sequence
-- subj_end - End coordinate of alignment on subject sequence
-- evalue - Expectation value
-- bit_score - Bit score of alignment
-- kog_id - KOG identifier
-- kog_name - KOG name
-- kog_length - Length of KOG consensus sequence

<taxon_oid>.pfam.tab.txt (from EBI's pfam_scan which uses HMMER 3.0)
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- query_start - Start coordinate of alignment on query gene
-- query_end - End coordinate of alignment on query gene
-- subj_start - Start coordinate of alignment on subject sequence
-- subj_end - End coordinate of alignment on subject sequence
-- evalue - Expectation value
-- bit_score - Bit score of alignment
-- pfam_id - Pfam identifier
-- pfam_name - Pfam name
-- pfam_length - Length of Pfam consensus sequence

<taxon_oid>.tigrfam.tab.txt (from hmmscan HMMER3.0)
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- query_start - Start coordinate of alignment on query gene
-- query_end - End coordinate of alignment on query gene
-- evalue - Expectation value
-- bit_score - Bit score of alignment
-- tigrfam_id - TIGRFAM identifier
-- tigrfam_name - TIGRFAM name

<taxon_oid>.ipr.tab.txt 
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- query_start - Start coordinate of alignment on query gene
-- query_end - End coordinate of alignment on query gene
-- domaindb - Original domain database.
-- domainid - ID on original domain database.
-- iprid - InterPro ID.
-- iprdesc - InterPro description.
-- go_info - Gene Ontology Information

<taxon_oid>.ko.tab.txt (from NCBI BLASTP on KEGG genes)
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- percent_identity - Perceent identity of aligned amino acid residues
-- query_start - Start coordinate of alignment on query gene
-- query_end - End coordinate of alignment on query gene
-- subj_start - Start coordinate of alignment on subject sequence
-- subj_end - End coordinate of alignment on subject sequence
-- evalue - Expectation value
-- bit_score - Bit score of alignment
-- ko_id - KEGG Orthology (KO) identifier
-- ko_name - KO name
-- EC - Enzyme Commission (EC) assignment from KO
-- img_ko_flag - 'Yes' (assigned by IMG pipeline); 'No' - from KEGG.

<taxon_oid>.signalp.tab.txt (SignalP)
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- feature_type - "cleavage"
-- start_coord - start coordinate of feature
-- end_coord - end coordinate of feature

<taxon_oid>.tmhmm.tab.txt (TMHMM)
-- gene_oid - Gene object identifier of query gene
-- gene_length - Length of protein sequence
-- feature_type - feature
-- start_coord - start coordinate of feature
-- end_coord - end coordinate of feature

<taxon_oid>.xref.tab.txt
-- gene_oid - Gene object identifier of query gene
-- db_name - External database
-- id - External ID corresponding to database

<taxon_oid>.crispr.txt (optional)
-- contig_id - Contig/Scaffold ID
-- crispr_no - CRISPR number
-- pos - Starting position of array element
-- repeat_seq - Repeat sequence
-- spacer_seq - Spacer sequence
-- tool_code - Single letter code for tool used

