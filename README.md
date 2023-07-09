# llama-magic


*****

This repository is outdated.  Please refer to https://github.com/FenyoLab/Llama_Nanobodies for the latest Llama Magic code!

*****

LM workflow

To create a database for searching :

1) dnatoprot_longest_nr.pl: translates to protein sequences, outputting only the longest non redundant reading frame  
input: fasta file from output of e.g. MiSeq (note, if paired end reads must be merged using a separate tool first)  
output: protein sequence fasta file

2) digest_fasta.pl: runs in silico digestion (trypsin) of the protein sequences, necessary for running Tandem search due to the high similarity of nanobody sequences  
input: output from the dnatoprot_longest_nr.pl script  
output: peptides fasta file

To run the search/mapping of peptides:

1) tandem.exe - the database for searching is the output of step 2 above

2) parse_xtandem_llama.pl: parses Tandem output xml file into txt file containing Tandem results  
input: directory of the result file "output.xml" from Tandem run, expectation threshold  
output: txt file with extracted information used in step 3 below

3) map_peptides_to_proteins.pl: the main LM algorithm - maps peptides back to sequences and finds CDR regions, ranks candidate nano bodies based on coverage of CDR regions, etc.  
input: file from step 2 above, fasta file containing protein sequences (longest_nr.fasta), xml file from Tandem output, 1/0 parameter whether to show the 'score' for a candidate nano body  
output: html file containing results of LM algorithm - the candidate nano bodies ranked * Note:  the matching peptides are shown with a link to a program called peptides.pl, which requires a web server to execute (it is a CGI script) and will show the details of the 'best matching spectrum' for that peptide.  
