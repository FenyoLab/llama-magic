#!/usr/local/bin/perl

#    map_peptides_to_proteins.pl
#
#    Copyright (C) 2014  Sarah Keegan
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
use strict;
#use warnings;
use POSIX;

my $img_source_plus = '/llama-magic-html/plus.gif';
my $img_source_minus = '/llama-magic-html/minus.gif';
my $img_source_star = '/llama-magic-html/greyplus.gif';

my $filename="";
my $fileroot;
my $filedir;
my $db_filename="";
my $tandem_output_filename="";
my $PEPTIDE_FILTER_MIN = 1.2; 
my $MIN_SEQ_LENGTH = 107; 
my $MIN_CDR3_COV_PERC = 25; #25; 
my $MAX_PROTEINS_IN_FILE = 10000;
my $MAX_GROUP_PROTEINS_IN_FILE = 10;
my $RK_FIXATION = 1;

my @nanobody_group; #the nanobodys are grouped by CDR3 region, allowing one mismatch

sub log10
{
	my $n = shift;
	return log($n)/log(10);
}

##
#changes for llama_magic web app:
#changed protein database to file name not directory (ARGV[1])
#put log file and output file in directory where tandem results txt file is located (ARGV[0]) rather than in current directory of script
#output html file name always includes a '1' even if there is only one output file (see function open_OUT_ALL)
#put link at bottom of output html file to go to the next output file (if it exists)
##
my $output_stats = 0; #if true, the program will out put the stats used for ranking in a separate CSV file for each sequence
my $stand_alone = 0; #this var will determine how we examine the input file containing XTamdem results - is it directly from XTAndem - i.e. run
		     #through llama-magic OR was it manually processed and provided to us - in each case the format is different
my $show_score = 0;
#if ($ARGV[0]=~/\w/) { $filename="$ARGV[0]"; } else { $stand_alone = 1; $filename="C:\\NCDIR\\Llama\\results\\19\\26\\tandem\\results\\output.xml.peptide_list.0.1.txt"; }
#if ($ARGV[1]=~/\w/) { $db_filename="$ARGV[1]"; } else { $stand_alone = 1; $db_filename="C:\\NCDIR\\Llama\\results\\19\\protein\\longest_nr.fasta"; }
#if ($ARGV[2]=~/\w/) { $tandem_output_filename="$ARGV[2]"; } else { $stand_alone = 1; $tandem_output_filename="C:\\NCDIR\\Llama\\results\\19\\26\\tandem\\results\\output.xml"; }

if ($ARGV[0]=~/\w/) { $filename="$ARGV[0]"; } else { $stand_alone = 1; $filename="./output.xml.peptide_list.0.1.txt"; }
if ($ARGV[1]=~/\w/) { $db_filename="$ARGV[1]"; } else { $stand_alone = 1; $db_filename="./longest_nr.fasta"; }
if ($ARGV[2]=~/\w/) { $tandem_output_filename="$ARGV[2]"; } else { $stand_alone = 1; $tandem_output_filename="./output.xml"; }
if ($ARGV[3]=~/\w/) { $show_score=$ARGV[3]; } 

my $cluster = 0;
if ($cluster)
{
	$stand_alone = 0;
	$show_score = 1;
	$output_stats = 1;
	$MIN_CDR3_COV_PERC = 1;
}

$filedir = $filename;
$filedir =~ s/[\/\\][^\/\\]+$//;

$filename =~ /[\/\\]([^\/\\]+)$/;
$fileroot = $1;
$fileroot =~ s/(\.\w\w\w)$//;

#open log file...
if(!open(LOG,">$filedir/$fileroot.count_proteins.log.txt"))
{
	print "Could not open for writing: $filedir/$fileroot.count_proteins.log.txt ($!)\n";
	exit(1);
}

#autoflush log file, so we an check progress more easily
select(LOG);
$|++; # autoflush LOG
select(STDOUT);

###########Reading in peptide data################################################################
my %peptides=(); #the peptides, $peptides{seq} = 1
my $count_pep=0; #the number of peptides found
my $count_pep_unique=0; #the number of unique peptides found

print LOG "Reading in peptide data...\n";

#read in database search results - the peptides that were found, and the expectation values
my $line="";
if(!open(IN,"$filename"))
{
	print LOG "Error: Could not open '$filename' ($!)\n";
	close(LOG);
	exit(1);
}

$line=<IN>;
while($line=<IN>)
{
	$line =~ s/\r\n$//;
	$line =~ s/\n$//;
	#chomp($line);
	
	if ((!$stand_alone && $line=~/^([A-Za-z]+)\t([0-9\-\+edED\.]*)\t([^\t]+)\t([^\t]+)\t([^\t]*)/) || #use this if we are running through llama magic
	    ($stand_alone && $line=~/^([A-Za-z]+)\t([0-9\-\+edED\.]*)/)) #use this if we don't have protein_uid and domain_id, and spectrum fields
	{
		my $original_pep=$1;
		my $expect=$2;
		
		my $protein_uid;
		my $domain_id;
		my $spectrum;
		
		if(defined $3) { $protein_uid = $3; }
		else { $protein_uid = ""; }
		if(defined $4) { $domain_id = $4; }
		else { $domain_id = ""; }
		if(defined $5) { $spectrum = $5; }
		else { $spectrum = ""; }
	
		my $pep = $original_pep;
		$pep=~tr/L/I/;
		if (!(defined $peptides{"\U$pep"})) 
		{#if peptide found more than once, only the lowest expectation value 
		 #and the associated original sequence is stored
			$count_pep_unique++; 
			if($spectrum) { $peptides{"\U$pep"}=[$expect,$original_pep,"$spectrum ($expect)", 1, $protein_uid, $domain_id]; }
			else { $peptides{"\U$pep"}=[$expect,$original_pep,"", 1, $protein_uid, $domain_id]; }
			
		}
		else
		{
			my $prev_expect = $peptides{"\U$pep"}[0];
			if($prev_expect > $expect)
			{
				$peptides{"\U$pep"}[0] = $expect;
				$peptides{"\U$pep"}[4] = $protein_uid;
				$peptides{"\U$pep"}[5] = $domain_id;
			}
			if($spectrum)
			{
				if($peptides{"\U$pep"}[2]) { $peptides{"\U$pep"}[2] .= "<br><br>$spectrum ($expect)"; }
				else { $peptides{"\U$pep"}[2] = "$spectrum ($expect)"; }
			}
			
			$peptides{"\U$pep"}[3]++;
		}
		$count_pep++;
	} 
	else { if ($line=~/\w/) { print LOG qq!Error: $line\n!; } }	
}
close(IN);
print LOG qq!$filename\n$count_pep peptides\n!;

#filter input peptides by expectation score
#filter: -log(e) * peptides_count >= 1.2 (peptides_count = # times peptide found in the file, log(e) = lowest expect value for peptide)
foreach (keys %peptides)
{
	if(abs($peptides{$_}[0] * $peptides{$_}[3]) < $PEPTIDE_FILTER_MIN)
	{#remove peptide from list
		delete $peptides{$_};
	}
}

###########Reading in protein data###############################################################
my %proteins=(); #the proteins, $proteins{name} = seq
my %protein_gene_counts;
my $count_total_seq = 0;
my $all_length_total_seq = 0;

print LOG "Reading in protein database data...\n";

if (!open (IN,"$db_filename"))
{
	print LOG "Error: Could not open '$db_filename' ($!)\n";
	close(LOG);
	exit(1);
}

print LOG "$db_filename\n";

my $name="";
my $gene_count;
my $sequence="";
my $count_seq=0; my $all_length_count_seq=0; 
while ($line=<IN>)
{
	my $name_; my $gene_count_ = 0;
	$line =~ s/\r\n$//;
	$line =~ s/\n$//;
	#chomp($line); 
	if ($line=~/^>(\S+) \+ (\d+) other/)
	{
		$name_=$1;
		$gene_count_ = $2 + 1;
		
	}
	elsif($line=~/^>(\S+)/)
	{
		$name_=$1;
		$gene_count_ = 1;
	}
	if($gene_count_ > 0)
	{#its a name line
		if ($name=~/\w/ and $sequence=~/\w/)
		{
			if(filter_input($sequence))
			{#check length from first M is atleast 100 and trim sequence 
				$proteins{$name}=$sequence;
				$protein_gene_counts{$name} = $gene_count;
				$count_seq++;
				$all_length_count_seq++;
			}
			else { $all_length_count_seq++; }
		}
		$name=$name_;
		$gene_count = $gene_count_;
		$sequence="";
	}
	else
	{#its a sequence line
		$sequence.="\U$line";
	}
	#if($count_seq > 10000) { last; }
}	
if ($name=~/\w/ and $sequence=~/\w/)
{
	if(filter_input($sequence))
	{#check length from first M is atleast 100		
		$proteins{$name}=$sequence;
		$protein_gene_counts{$name} = $gene_count;
		$count_seq++;
		$all_length_count_seq++;
	}
	else { $all_length_count_seq++; }
}
$count_total_seq += $count_seq;
$all_length_total_seq += $all_length_count_seq;
close(IN);
print LOG qq!$count_seq sequences\n!;

print LOG "Found $count_total_seq out of $all_length_total_seq sequences with minimum length.\n";


################Find the 3 CDR regions and store their positions:#################################################
my %CDR1; #$CDR1{'prot-name'} = [pos_st_cdr, pos_end_cdr] (array ref)
my %CDR2; #$CDR2{'prot-name'} = [pos_st_cdr, pos_end_cdr] (array ref)
my %CDR3; #$CDR3{'prot-name'} = [pos_st_cdr, pos_end_cdr] (array ref)
my %CDR3_seq;
#
#print "Finding CDR regions...\n";
#foreach my $name (keys %proteins)
#{
#	my $seq = $proteins{$name};
#	my $l_pos; my $r_pos;
#	
#	find_cdr1($seq, $l_pos, $r_pos);
#	$CDR1{$name} = [$l_pos, $r_pos];
#	
#	find_cdr2($seq, $l_pos, $r_pos);
#	$CDR2{$name} = [$l_pos, $r_pos];
#	
#	find_cdr3($seq, $l_pos, $r_pos);
#	$CDR3{$name} = [$l_pos, $r_pos];
#	$CDR3_seq{$name} = substr($seq, $l_pos, $r_pos-$l_pos+1);
#}	

#################Match peptides to proteins#########################################################################
my %peptide_proteins; #the proteins found (in the fasta dbs) that contain this peptide - peptide_proteins{pep-seq} = #prot-name1##prot-name2#...
my %peptide_proteins_count; #the number of proteins found that match this peptide - protein_peptides_count{pep-seq} = n
my %protein_peptides; #the peptides found in this protein - protein_peptides{prot-name} - #pep-seq1##pep-seq2#...
my %protein_peptides_count; #the number of peptides found in this protein - protein_peptides_count{prot-name} = n
my %protein_peptides_pos;
my %protein_total_coverage_by_aa;
my %protein_cdr1_cover;
my %protein_cdr2_cover;
my %protein_cdr3_cover;
my %protein_total_cover;
my @protein_coverage_stats;

print LOG "Matching peptides to proteins...\n";

if ($output_stats)
{
	if(!open(STAT,">$filedir/$fileroot.count_proteins.stat.csv"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.count_proteins.stat.txt ($!)\n";
		$output_stats = 0;
	}
	else
	{
		print STAT "ID,CDR3_PERC,CDR3_LEN,HT_SEQ_COUNT,CDR123_PERC,CDR123_LEN,SEQ_PERC,SEQ_LEN,CDR1_PERC,CDR1_LEN,CDR2_PERC,CDR2_LEN,SEQ\n";
	}
}

my $count_proteins = 0; my $stat_i = 0;
#foreach my $pep (keys %peptides) { $peptide_proteins_count{$pep}=0; }
foreach my $name (keys %proteins)
{
	#$protein_peptides_count{$name} = 0;
	my $seq=$proteins{$name};
	$seq=~tr/L/I/;
	my %protein_coverage;
	foreach my $pep (keys %peptides)
	{
		my $pos = 0; my $num_matches = 0;
		while(($pos = index($seq, $pep, $pos)) >= 0)
		{
			my $char = '';
			if($pos >= 0) { $char = substr($seq, $pos-1, 1); }
			if(!$RK_FIXATION || ($char eq 'K' || $char eq 'R')) #not allowing peptides at beginning of sequence, MUST have R/K before peptide (simplification)
			{
				$num_matches++;
				if($num_matches > 1) { print LOG "Note: peptide found more than once in a protein: '$name', '$pep'\n"; }
				$protein_peptides{$name}.="#$pep#";
				$protein_peptides_count{$name}++;
				$protein_peptides_pos{$name}.="#$pos#";
				
				#collect 'uniqueness' information (all proteins for a peptide)
				$peptide_proteins{$pep}.="$name<br>";
				$peptide_proteins_count{$pep}++;
				
				#tally peptide coverage for this protein
				my $pep_len = length($pep);
				for(my $i = $pos; $i < ($pos+$pep_len); $i++) { $protein_coverage{$i} = 1; }
				
			}
			$pos = $pos + length($pep);
			my $temp;
		}
	}
	
	my $cdr1_l_pos; my $cdr1_r_pos;	my $cdr2_l_pos; my $cdr2_r_pos;	my $cdr3_l_pos; my $cdr3_r_pos;	
	find_cdr1($proteins{$name}, $cdr1_l_pos, $cdr1_r_pos);
	find_cdr2($proteins{$name}, $cdr2_l_pos, $cdr2_r_pos);
	find_cdr3($proteins{$name}, $cdr3_l_pos, $cdr3_r_pos);
	
	#calculate coverage percent for this protein (for sorting and display)
	my $cover_n = 0; my $cover_cdr = 0; my $cover_cdr3 = 0; my $cover_cdr2 = 0; my $cover_cdr1 = 0;
	foreach (keys %protein_coverage)
	{
		if($cdr1_l_pos <= $_ && $cdr1_r_pos >= $_) { $cover_cdr++; $cover_cdr1++; }
		if($cdr2_l_pos <= $_ && $cdr2_r_pos >= $_) { $cover_cdr++; $cover_cdr2++; }
		if($cdr3_l_pos <= $_ && $cdr3_r_pos >= $_) { $cover_cdr++; $cover_cdr3++; }
		$cover_n++;
	}
	#total 3 cdrs length
	my $cdr_len = ($cdr1_r_pos - $cdr1_l_pos + 1) + ($cdr2_r_pos - $cdr2_l_pos + 1) + ($cdr3_r_pos - $cdr3_l_pos + 1);
	my $seq_len = length($seq);
	my $cdr3_len = $cdr3_r_pos - $cdr3_l_pos + 1;
	my $cdr3_cover_perc = ($cover_cdr3 / $cdr3_len) * 100;
	#if(($MIN_CDR3_COV_PERC == 0 && $cdr3_cover_perc > $MIN_CDR3_COV_PERC) ||
	#   ($MIN_CDR3_COV_PERC > 0 && $cdr3_cover_perc >= $MIN_CDR3_COV_PERC))
	if ($cdr3_cover_perc >= $MIN_CDR3_COV_PERC) 
	{
		#create the sorting array:
		$protein_coverage_stats[$stat_i][0] = sprintf("%.1f", $cdr3_cover_perc);
		$protein_coverage_stats[$stat_i][1] = $cdr3_len;
		$protein_coverage_stats[$stat_i][2] = $protein_gene_counts{$name};
		$protein_coverage_stats[$stat_i][3] = sprintf("%.1f", ($cover_cdr / $cdr_len) * 100);
		$protein_coverage_stats[$stat_i][4] = $cdr_len;
		$protein_coverage_stats[$stat_i][5] = sprintf("%.1f", ($cover_n / $seq_len) * 100);
		$protein_coverage_stats[$stat_i][6] = $seq_len;
		$protein_coverage_stats[$stat_i][7] = $name;
		my $s = 8*($cover_cdr3 / $cdr3_len) +
			2*($cover_cdr2 / ($cdr2_r_pos - $cdr2_l_pos + 1)) +
			2*($cover_cdr1 / ($cdr1_r_pos - $cdr1_l_pos + 1)) +
			2*($cover_n / $seq_len) +
			$cdr3_len/15 +
			($cdr2_r_pos - $cdr2_l_pos + 1)/10 +
			($cdr1_r_pos - $cdr1_l_pos + 1)/9;
		$protein_coverage_stats[$stat_i][8] = sprintf("%.1f", $s);
		$stat_i++;
		
		#store more info about the protein for display:
		$protein_total_coverage_by_aa{$name} = \%protein_coverage;
		$protein_cdr1_cover{$name} = $cover_cdr1;
		$protein_cdr2_cover{$name} = $cover_cdr2;
		$protein_cdr3_cover{$name} = $cover_cdr3;
		$protein_total_cover{$name} = $cover_n;
		
		$CDR1{$name} = [$cdr1_l_pos, $cdr1_r_pos];
		$CDR2{$name} = [$cdr2_l_pos, $cdr2_r_pos];
		$CDR3{$name} = [$cdr3_l_pos, $cdr3_r_pos];
		$CDR3_seq{$name} = substr($seq, $cdr3_l_pos, $cdr3_r_pos-$cdr3_l_pos+1);
	}
	
	if ($output_stats)
	{
		my $cdr1_len = $cdr1_r_pos - $cdr1_l_pos + 1;
		my $cdr2_len = $cdr2_r_pos - $cdr2_l_pos + 1;
		my $cdr1_perc = sprintf("%.1f", ($cover_cdr1 / $cdr1_len) * 100);
		my $cdr2_perc = sprintf("%.1f", ($cover_cdr2 / $cdr2_len) * 100);
				
		print STAT "$name," .
			sprintf("%.1f", $cdr3_cover_perc) .
			",$cdr3_len,$protein_gene_counts{$name}," .
			sprintf("%.1f", ($cover_cdr / $cdr_len) * 100) .
			",$cdr_len," .
			sprintf("%.1f", ($cover_n / $seq_len) * 100) .
			",$seq_len,$cdr1_perc,$cdr1_len,$cdr2_perc,$cdr2_len,$proteins{$name}\n";
	}
	
	if($count_proteins % 50000 == 0) { print LOG qq!$count_proteins ($count_total_seq)\n!; }
	$count_proteins++;
}

if ($output_stats) { close(STAT); }

foreach my $pep (keys %peptides) { $peptide_proteins{$pep} =~ s/, $//; } #take off commas at end of name list

print LOG "Done!\n";


############Sort proteins based on CDR coverage, length, etc.################################################################
print LOG "Sorting the proteins for ouput...\n";

my @proteins_to_sort;
my $proteins_to_sort_count=0;
my @proteins_sorted;

if ($show_score)
{
	@proteins_sorted = sort { $b->[8] <=> $a->[8] || $b->[2] <=> $a->[2] } @protein_coverage_stats;
}
else
{
	@proteins_sorted = sort
	{ $b->[0] <=> $a->[0] || $b->[1] <=> $a->[1] || $b->[2] <=> $a->[2] || $b->[3] <=> $a->[3] || $b->[4] <=> $a->[4] || $b->[5] <=> $a->[5] || $b->[6] <=> $a->[6] }
	@protein_coverage_stats;
}

$proteins_to_sort_count = $#proteins_sorted + 1;

print LOG "Proteins sorted! ($proteins_to_sort_count)\n";

############Group the proteins by CDR3 region w/ one AA difference allowed, combine groups###################################
print LOG "Grouping proteins on CDR3 region...\n";

# SLOW WAY ##########
#print "Grouping proteins on CDR3 region...\n";
#
##group the proteins based on the CDR3 region - the same or 1 AA difference is allowed
##make the groups:
#for(my $i = 0; $i <= $#proteins_sorted; $i++)
#{
#	my $cdr3_seq_i = $CDR3_seq{$proteins_sorted[$i][6]};
#	for(my $j = $i+1; $j <= $#proteins_sorted; $j++)
#	{
#		my $cdr3_seq_j = $CDR3_seq{$proteins_sorted[$j][6]};
#		if(compare_cdrs($cdr3_seq_i, $cdr3_seq_j))
#		{
#			${$nanobody_group[$i]}{$j} = 1;
#		}
#	}
#	if(!$nanobody_group[$i]) { $nanobody_group[$i] = (); }
#	if($i % 50 == 0) { print "$i...\n"; }
#}
#
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	print LOG "$i: ";
#	if(defined ${$nanobody_group[$i]}{-1}) { print LOG "x"; }
#	else { print LOG join ', ', sort keys %{$nanobody_group[$i]}; }
#	print LOG "\n";
#}
#
#print "Consolidating groups...\n"; 
#
##consolidate the groups - if 2 or more groups share the same sequence, combine them
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	if($nanobody_group[$i])
#	{
#		for(my $j = $i+1; $j <= $#nanobody_group; $j++)
#		{
#			consolidate_groups($i, $j);
#			print "$j...\n";
#		}
#	}
#	#if($i % 50 == 0) { print "$i...\n"; }
#	print "$i...\n";
#}
#
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	print LOG "$i: ";
#	if(defined ${$nanobody_group[$i]}{-1}) { print LOG "x"; }
#	else { print LOG join ', ', sort keys %{$nanobody_group[$i]}; }
#	print LOG "\n";
#}

my %lowest_matching_i;
for(my $i = 0; $i <= $#proteins_sorted; $i++)
{
	if($i % 1000 == 0) { print LOG "$i...\n"; }
	if(defined ${$nanobody_group[$i]}{-1}) { next; }
	my %matching_is; my @matching_js;
	my $cdr3_seq_i = $CDR3_seq{$proteins_sorted[$i][7]};
	for(my $j = $i+1; $j <= $#proteins_sorted; $j++)
	{
		if(defined ${$nanobody_group[$j]}{-1}) { next; }
		my $cdr3_seq_j = $CDR3_seq{$proteins_sorted[$j][7]};
		if(compare_cdrs($cdr3_seq_i, $cdr3_seq_j))
		{#if cdrs match...
			if(defined $lowest_matching_i{$j})
			{#if protein j already matches a previous protein, record this protein's position 
				$matching_is{$lowest_matching_i{$j}} = 1;
			}
			push @matching_js, $j; #add protein j to list of matches to protein i
		}
	}
	if($#matching_js >= 0)
	{#if there are proteins that matched to protein i
		if(keys %matching_is)
		{#if there are other proteins that matched to the matches of protein i
			
			#get lowest group and add matching js
			my @sorted_matching_is = sort {$a <=> $b} keys %matching_is;
			my $lowest_i = $sorted_matching_is[0]; #the lowest matching protein group is the group to put all matches into
			
			foreach my $cur_j (@matching_js)
			{#for all matching j's to protein i, put them in the lowest i group
				${$nanobody_group[$lowest_i]}{$cur_j} = 1; #mark j as being in lowest matching protein group
				$lowest_matching_i{$cur_j} = $lowest_i; #set lowest matching i for this j
				${$nanobody_group[$cur_j]}{-1} = 1; #protein j is in another group, x out its entry
			}
			
			foreach my $cur_i (@sorted_matching_is)
			{#for each of the other proteins found that matched proteins matching protein i, we need to put them in the lowest i group
				if($cur_i != $lowest_i)
				{
					foreach my $cur_j (keys %{$nanobody_group[$cur_i]})
					{
						${$nanobody_group[$lowest_i]}{$cur_j} = 1;
						$lowest_matching_i{$cur_j} = $lowest_i;
					}
					${$nanobody_group[$cur_i]}{-1} = 1;
				}
			}
			
			#set the group for protein i as well
			${$nanobody_group[$lowest_i]}{$i} = 1;
			$lowest_matching_i{$i} = $lowest_i;
		}
		else
		{#its just new proteins that match, no groups
			foreach my $cur_j (@matching_js)
			{
				${$nanobody_group[$i]}{$cur_j} = 1; #mark protein j as being in protein i's group
				$lowest_matching_i{$cur_j} = $i; #set lowest matching i to the current i for protein j
				${$nanobody_group[$cur_j]}{-1} = 1; #protein j is in another group, x out its entry
			}
		}
	}
	else
	{
		$nanobody_group[$i] = ();
	}
}

#printing Stats and verifying groups...
my $num = $#nanobody_group+1;
print LOG "nanobody_group has $num members.\n";
$num = $#proteins_sorted+1;
print LOG "proteins_sorted has $num members.\n";
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	print LOG "$i: ";
#	if(defined ${$nanobody_group[$i]}{-1}) { print LOG "x"; }
#	else { print LOG join ', ', sort { $a <=> $b } keys %{$nanobody_group[$i]}; }
#	print LOG "\n";
#}

#do some checking
#print "Verifying groups...\n";
#my @checking_group;
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	if(!${$nanobody_group[$i]}{-1})
#	{
#		$checking_group[$i]++;
#		foreach (keys %{$nanobody_group[$i]})
#		{
#			$checking_group[$_]++;
#		}
#		
#	}
#}
#for(my $i = 0; $i <= $#nanobody_group; $i++)
#{
#	if($checking_group[$i] != 1) { print "position $i of checking_group is checking_group[$i].\n"}
#}

###### Outputting results to a txt file ############### <- only outputs data for > min coverage of CDR3 region - output all data above

if ($output_stats)
{
	#open log file...
	if(!open(STAT,">$filedir/$fileroot.count_proteins.group-stat.csv"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.count_proteins.stat.txt ($!)\n";
	}
	else
	{
		print STAT "ID,GROUP_HEAD,CDR3_PERC,CDR3_LEN,HT_SEQ_COUNT,CDR123_PERC,CDR123_LEN,SEQ_PERC,SEQ_LEN,CDR1_PERC,CDR1_LEN,CDR2_PERC,CDR2_LEN,SEQ\n";
		for(my $i = 0; $i <= $#nanobody_group; $i++)
		{
			
			if(! (defined ${$nanobody_group[$i]}{-1})) 
			{
				my $name = $proteins_sorted[$i][7];
				#print out stats for head of group
				
				###
				my $cdr1_len = $CDR1{$name}[1] - $CDR1{$name}[0] + 1; my $cdr2_len = $CDR2{$name}[1] - $CDR2{$name}[0] + 1;
				my $cdr1_perc = sprintf("%.1f", ($protein_cdr1_cover{$name} / $cdr1_len) * 100);
				my $cdr2_perc = sprintf("%.1f", ($protein_cdr2_cover{$name} / $cdr2_len) * 100);
				print STAT "$name,1,$proteins_sorted[$i][0],$proteins_sorted[$i][1],$proteins_sorted[$i][2],$proteins_sorted[$i][3],$proteins_sorted[$i][4],$proteins_sorted[$i][5],$proteins_sorted[$i][6],$cdr1_perc,$cdr1_len,$cdr2_perc,$cdr2_len\n";
				
				#print out stats for group members
				my @current_group = sort { $a <=> $b } keys %{$nanobody_group[$i]};
				for(my $k = 0; $k <= $#current_group; $k++)
				{
					my $j = $current_group[$k];
					$name = $proteins_sorted[$j][7];
					
					###
					$cdr1_len = $CDR1{$name}[1] - $CDR1{$name}[0] + 1; my $cdr2_len = $CDR2{$name}[1] - $CDR2{$name}[0] + 1;
					$cdr1_perc = sprintf("%.1f", ($protein_cdr1_cover{$name} / $cdr1_len) * 100);
					$cdr2_perc = sprintf("%.1f", ($protein_cdr2_cover{$name} / $cdr2_len) * 100);
					print STAT "$name,0,$proteins_sorted[$j][0],$proteins_sorted[$j][1],$proteins_sorted[$j][2],$proteins_sorted[$j][3],$proteins_sorted[$j][4],$proteins_sorted[$j][5],$proteins_sorted[$i][6],$cdr1_perc,$cdr1_len,$cdr2_perc,$cdr2_len,$proteins{$name}\n";
				
				}
				
			}
		}
		close(STAT);
	}
}

#return 0;

##############Output results to the html file###################################################################################

my %protein_peptides_unique_count;
my %peptides_done;
my $proteins_unique_count=0;

print LOG "Outputting results to the file ($proteins_to_sort_count proteins with > $MIN_CDR3_COV_PERC% CDR3 coverage)...\n";

open_OUT_ALL(1);

if($proteins_to_sort_count == 0)
{
	print OUT_ALL "No matches found.\n";
	close_OUT_ALL(0, 1); 
	close(LOG);
	exit(0);
}


my $group_div = 0; my @current_group; my $group_div_head; my $indent = ""; my $make_new_file = 0; my $file_num = 1; my $num_group_members;
for(my $group_count = -1, my $cur_group_i = 0, my $real_group_count = 1, my $total_proteins_in_file = 0; $group_count <= $#nanobody_group; $cur_group_i++, $total_proteins_in_file++)
{
	my $proteins_sorted_count; my $group_head;
	
	if($total_proteins_in_file > $MAX_PROTEINS_IN_FILE)
	{ $total_proteins_in_file = 0; $make_new_file = 1; }
	
	if($cur_group_i > $#current_group || ($cur_group_i+1) >= $MAX_GROUP_PROTEINS_IN_FILE) #$cur_group_i starts at -1
	{#go to the next group
		if($cur_group_i <= $#current_group) { print OUT_ALL '          ...<br>'; } #there's more in the group, show ...
		$group_count++;
		
		if($group_count > $#nanobody_group) { last; } #this isn't a new group, it's the end!
		#if($group_count > 30) { last; } 
		
		if(defined ${$nanobody_group[$group_count]}{-1})
		{#next group doesn't exist (it's part of another group), skip
			@current_group = (); $cur_group_i = 0;
			$total_proteins_in_file--;
			next;
		}
		else
		{
			$proteins_sorted_count = $group_count;
			@current_group = sort { $a <=> $b } keys %{$nanobody_group[$group_count]};
			$num_group_members = $#current_group+2; #add in an extra 1 for the group head
			$group_head = 1; $cur_group_i = -1;
		}
		
		if($make_new_file)
		{
			close_OUT_ALL($group_div, ++$file_num);
			open_OUT_ALL($file_num);
			$make_new_file = 0;
		}
	}
	else
	{
		$proteins_sorted_count = $current_group[$cur_group_i];
		$group_head = 0;
	}
	
	#read in the name of the next protein in the sorted array
	my @protein_peptides; 
	my $num_pp = 0;
	my $name;
	$name = $proteins_sorted[$proteins_sorted_count][7];
	$protein_peptides_unique_count{$name}=0;
	my $temp=$protein_peptides{$name};
	my $temp1 = $protein_peptides_pos{$name};
	while($temp=~s/^#([^#]+)#//)
	{#read in peptides for this protein
		my $pep=$1;
		$protein_peptides[$num_pp][0] = $pep;
		$temp1 =~ s/^#([^#]+)#//;
		$protein_peptides[$num_pp][1] = $1; #the starting position of the peptide in the protein
		$num_pp++
	}
	#sort peptides based on starting position in the protein
	my @sorted_protein_peptides = sort { $a->[1] <=> $b->[1] } @protein_peptides;
	 
	#print out the protein with stats
	#format: CDR1: X% (X/X); CDR2: X% (X/X); CDR3: X% (X/X); combined CDR: X% (X/X); overall: X% (X/X)
	my $cdr1_len = $CDR1{$name}[1] - $CDR1{$name}[0] + 1; my $cdr2_len = $CDR2{$name}[1] - $CDR2{$name}[0] + 1;
	my $cdr1_perc = sprintf("%.1f", ($protein_cdr1_cover{$name} / $cdr1_len) * 100);
	my $cdr2_perc = sprintf("%.1f", ($protein_cdr2_cover{$name} / $cdr2_len) * 100);
	my $total_cdr_cover = $protein_cdr1_cover{$name} + $protein_cdr2_cover{$name} + $protein_cdr3_cover{$name};
	#my $score = 8*($protein_cdr3_cover{$name}/$proteins_sorted[$proteins_sorted_count][1]) + 2*($protein_cdr2_cover{$name} / $cdr2_len) +
	#	     2*($protein_cdr1_cover{$name} / $cdr1_len) + 2*($protein_total_cover{$name}/$proteins_sorted[$proteins_sorted_count][6]) + 
	#	     $proteins_sorted[$proteins_sorted_count][1]/15 + $cdr2_len/10 + $cdr1_len/9;
	#my $score = sprintf("%.2f", $proteins_sorted[$proteins_sorted_count][8]);
	if($group_head)
	{
		if($group_div) { print OUT_ALL qq!</div>!; $group_div = 0; }
		print OUT_ALL "<br><b title='rank of the group'>$real_group_count</b>";
		
		#add spaces after rank so all results line up
		my $num_sp = 3 - floor(log10($real_group_count)+1);
		for(my $s = 0; $s < $num_sp; $s++) { print OUT_ALL ' ';  }
	
		$real_group_count++;
		if($#current_group >= 0)
		{
			$group_div_head = $proteins_sorted_count;	
			print OUT_ALL qq!<img src="$img_source_plus" title="click to expand the group for viewing more sequences" onclick="ec('d_g_$proteins_sorted_count', 'i_g_$proteins_sorted_count')" id="i_g_$proteins_sorted_count" style="cursor:hand;" alt="+" />!;
		}
		else
		{
			print OUT_ALL qq!<img src="$img_source_star" style="cursor:hand;" alt="*" disabled="true" title="no more sequences in this group"/>!;
		}
	}
	elsif($cur_group_i == 0)
	{#its the first sequence in group, print out div before it
		print OUT_ALL qq!<div id="d_g_$group_div_head" style="display:none">!;
		$group_div = 1;
	}
	
	if ($group_head)
	{
		if ($num_group_members > $MAX_GROUP_PROTEINS_IN_FILE) { print OUT_ALL "<b title='number of sequences in this group (first $MAX_GROUP_PROTEINS_IN_FILE are shown)'>($num_group_members)</b>"; }
		else { print OUT_ALL "<b title='number of sequences in this group'>($num_group_members)</b>"; }
		
		#add spaces after num group members so all results line up
		my $num_sp = 3 - floor(log10($num_group_members)+1);
		for(my $s = 0; $s < $num_sp; $s++) { print OUT_ALL ' '; }
	}
	else
	{#print out indentation so group members line up with head of group
		print OUT_ALL '          ';
	}
	print OUT_ALL qq!<img src="$img_source_plus" title="click to expand the sequence for peptide mapping information" onclick="ec('d_$proteins_sorted_count', 'i_$proteins_sorted_count')" id="i_$proteins_sorted_count" style="cursor:hand;" alt="+" />!;
	if ($show_score)
	{
		print OUT_ALL qq!<b title="sequence name">$name</b> <b title="SCORE=8*CDR3-COV + 2*CDR2-COV + 2*CDR1-COV + 2*SEQ-COV + CDR3-LEN/15 + CDR2-LEN/10 + CDR1-LEN/9">Score: $proteins_sorted[$proteins_sorted_count][8];</b> <b title="CDR1 coverage">CDR1: $cdr1_perc\% ($protein_cdr1_cover{$name}/$cdr1_len);</b> <b title="CDR2 coverage">CDR2: $cdr2_perc\% ($protein_cdr2_cover{$name}/$cdr2_len);</b> <b title="CDR3 coverage">CDR3: $proteins_sorted[$proteins_sorted_count][0]\% ($protein_cdr3_cover{$name}/$proteins_sorted[$proteins_sorted_count][1]);</b> <b title="CDR 1, 2 and 3 combined coverage">combined CDR: $proteins_sorted[$proteins_sorted_count][3]\% ($total_cdr_cover/$proteins_sorted[$proteins_sorted_count][4]);</b> <b title="overall sequence coverage">overall: $proteins_sorted[$proteins_sorted_count][5]\% ($protein_total_cover{$name}/$proteins_sorted[$proteins_sorted_count][6]);</b> <b title="number of HT-DNA sequencing reads that produced this sequence">HT-seq count: $proteins_sorted[$proteins_sorted_count][2]</b>\n             !;
	}
	else
	{
		print OUT_ALL qq!<b title="sequence name">$name</b> <b title="CDR1 coverage">CDR1: $cdr1_perc\% ($protein_cdr1_cover{$name}/$cdr1_len);</b> <b title="CDR2 coverage">CDR2: $cdr2_perc\% ($protein_cdr2_cover{$name}/$cdr2_len);</b> <b title="CDR3 coverage">CDR3: $proteins_sorted[$proteins_sorted_count][0]\% ($protein_cdr3_cover{$name}/$proteins_sorted[$proteins_sorted_count][1]);</b> <b title="CDR 1, 2 and 3 combined coverage">combined CDR: $proteins_sorted[$proteins_sorted_count][3]\% ($total_cdr_cover/$proteins_sorted[$proteins_sorted_count][4]);</b> <b title="overall sequence coverage">overall: $proteins_sorted[$proteins_sorted_count][5]\% ($protein_total_cover{$name}/$proteins_sorted[$proteins_sorted_count][6]);</b> <b title="number of HT-DNA sequencing reads that produced this sequence">HT-seq count: $proteins_sorted[$proteins_sorted_count][2]</b>\n             !;
	}
	
	#underline AA's that correspond to identified peptides:
	#show the 3 CDR regions in red
	my $cov_by_aa_ref = $protein_total_coverage_by_aa{$name};
	my @seq = split('', $proteins{$name});
	my $u_on = 0;
	my $cdr_i = 1;
	for(my $i = 0; $i <= $#seq; $i++)
	{
		my $pre = ""; my $post = "";
		
		if(($i == $CDR1{$name}[0]) || ($i == $CDR2{$name}[0]) || ($i == $CDR3{$name}[0]))
		{ $pre .= "<font title=\"CDR $cdr_i\" color=\"\#FF0000\">"; $cdr_i++;}
		if(($i == $CDR1{$name}[1]) || ($i == $CDR2{$name}[1]) || ($i == $CDR3{$name}[1]))
		{ $post .= "</font>"; }
		
		if(defined $$cov_by_aa_ref{$i}) { if(!$u_on) { $pre .= "<u title='peptide coverage'>"; $u_on = 1; } }
		elsif($u_on) { $pre .= "</u>"; $u_on = 0; }
		
		print OUT_ALL "$pre$seq[$i]$post";	
	}
	if($u_on) { print OUT_ALL "</u>"; }	
	print OUT_ALL "\n"; 

	#print out all sorted peptides
	print OUT_ALL qq!<div id="d_$proteins_sorted_count" style="display:none">!; 
	for(my $i = 0; $i <= $#sorted_protein_peptides; $i++) 
	{
		my $expect = $peptides{$sorted_protein_peptides[$i][0]}[0];
		my $orig = $peptides{$sorted_protein_peptides[$i][0]}[1];
		
		if($i > 0) { print OUT_ALL "\n"; }
		printf OUT_ALL "   <b title='the best score (log(e)) for this peptide'>%-9s</b> ", $expect;
		
		for(my $s = 0; $s < $sorted_protein_peptides[$i][1]; $s++) { print OUT_ALL " "; }
		
		my $source = $peptides{$sorted_protein_peptides[$i][0]}[2];
		my $num_found = $peptides{$sorted_protein_peptides[$i][0]}[3];
		my $protein_uid = $peptides{$sorted_protein_peptides[$i][0]}[4];
		my $domain_id = $peptides{$sorted_protein_peptides[$i][0]}[5];
		my $num_proteins = $peptide_proteins_count{$sorted_protein_peptides[$i][0]};
		my $protein_names = $peptide_proteins{$sorted_protein_peptides[$i][0]};
	
		#score:
		my $score = sprintf("%.1f", $num_proteins/$num_group_members);
		print OUT_ALL qq!<a title="click here to view best matching spectrum" href="/llama-magic-cgi/peptide.pl?path=$tandem_output_filename&uid=$protein_uid&id=$domain_id&label=$orig">$orig</a>!;
		print OUT_ALL " <b title='uniqueness score (# matching sequences/# group members)'>($score)</b>  <b title='spectra count'>($num_found)</b>";
		print OUT_ALL qq!<img src="$img_source_plus" title="click to expand for spectra information" onclick="ec('d_$proteins_sorted_count-$i', 'i_$proteins_sorted_count-$i')" id="i_$proteins_sorted_count-$i" style="cursor:hand;" alt="+" />!;
		print OUT_ALL qq!<div id="d_$proteins_sorted_count-$i" style="display:none">$source</div>!;
		if($num_proteins < 100)
		{
			print OUT_ALL "  <b title='number of sequences matched for this peptide across entire database'>($num_proteins)</b>";
			print OUT_ALL qq!<img src="$img_source_plus" title="click to expand for all matched sequences" onclick="ec('d_$proteins_sorted_count-$i-2', 'i_$proteins_sorted_count-$i-2')" id="i_$proteins_sorted_count-$i-2" style="cursor:hand;" alt="+" />!;
			print OUT_ALL qq!<div id="d_$proteins_sorted_count-$i-2" style="display:none">$protein_names</div>!;
		}
		else
		{
			print OUT_ALL "  <b title='number of sequences matched for this peptide across entire database'>($num_proteins)</b>";
			print OUT_ALL qq!<img src="$img_source_star" style="cursor:hand;" alt="*" disabled="true" title="100 or more sequences matched this peptide"/>!;
		}	
	}
	print OUT_ALL qq!</div>!; 
	
}

close_OUT_ALL($group_div, 1); 

close(LOG);

##############subroutines##########################################################################################
sub filter_input
{#trim until the first M and check if length is >= 107
 #also check if 'QVT' is is at the end (last 10 AA's of sequence)
	my $seq = $_[0];
	my $m_i = index($seq, 'M');
	if($m_i >= 0) { $_[0] = substr($seq, $m_i); } else { $_[0] = ""; }
	
	if(length($_[0]) >= $MIN_SEQ_LENGTH)
	{
		#length okay, 'M' found, check 'QVT'
		my $end_seq = substr($_[0], length($_[0])-10); #look n-10 to n-1 (0-based)
		if($end_seq =~ /QVT/) { return 1; }
		return 0;
	}
	
	return 0; 
}

sub find_cdr1
{#finds the cdr1 region and sets args 2, 3 to starting and ending pos.
	my $seq = $_[0];
	
	#STARTING POS. OF CDR1:
	my $left_area = substr($seq, 20, 7); #look from pos. 20 - 26 of seq (0-based)
	my $la_i; my $left_cdr = -1;
	if(($la_i = index($left_area, 'SC')) < 0)
	{#didn't find 'SC', look for 'C'
		$la_i = index($left_area, 'C');
	} else { $la_i++; } #'C' is our marker, so advance past 'S'
	if($la_i >= 0) { $left_cdr = $la_i + 20 + 5; } #CDR1 starts at 'C' + 5 (add 20 to put it back in the full sequence)
	
	#ENDING POS. OF CDR1:
	my $right_area = substr($seq, 32, 9); #look from pos. 32 - 40 of seq (0-based)
	my $ra_i; my $right_cdr = -1;
	if($right_area =~ /(W.R)/)
	{#if we found 'WXR', find its index
		$ra_i = index($right_area, $1);
	}
	else { $ra_i = index($right_area, 'W'); } #didn't find 'WXR', look for 'W'
	if($ra_i >= 0) { $right_cdr = $ra_i + 32 - 1; } #CDR1 ends at 'W' - 1 (add 32 to put it back in the full sequence)
	
	#check if st/end found and if not follow rules:
	if($left_cdr == -1 && $right_cdr == -1) { $left_cdr = 28; $right_cdr = 36; }
	elsif($left_cdr == -1) { $left_cdr = $right_cdr - 8; }
	elsif($right_cdr == -1) { $right_cdr = $left_cdr + 8; }
	
	$_[1] = $left_cdr;
	$_[2] = $right_cdr;
	
	return 1;
}

sub find_cdr2
{#finds the cdr2 region and sets args 2, 3 to starting and ending pos.
	my $seq = $_[0];
	
	#STARTING POS. OF CDR2:
	my $left_area = substr($seq, 32, 9); #look from pos. 32 - 40 of seq (0-based)
	my $la_i; my $left_cdr = -1;
	if($left_area =~ /(W.R)/)
	{#if we found 'WXR', find its index
		$la_i = index($left_area, $1);
	}
	else { $la_i = index($left_area, 'W'); } #didn't find 'WXR', look for 'W'
	if($la_i >= 0) { $left_cdr = $la_i + 32 + 14; } #CDR2 starts at 'W' + 14 (add 32 to put it back in the full sequence)
	
	#ENDING POS. OF CDR2:
	my $right_area = substr($seq, 63, 10); #look from pos. 63 - 72 of seq (0-based)  
	my $ra_i; my $right_cdr = -1;
	$ra_i = index($right_area, 'RF');
	if($ra_i >= 0) { $right_cdr = $ra_i + 63 - 8; } #CDR2 ends at 'R' - 8 (add 63 to put it back in the full sequence)
	
	#check if st/end found and if not follow rules:
	if($left_cdr == -1 && $right_cdr == -1) { $left_cdr = 51; $right_cdr = 60; }
	elsif($left_cdr == -1) { $left_cdr = $right_cdr - 9; }
	elsif($right_cdr == -1) { $right_cdr = $left_cdr + 9; }
	
	$_[1] = $left_cdr;
	$_[2] = $right_cdr;
	
	return 1;
}

sub find_cdr3
{#finds the cdr3 region and sets args 2, 3 to starting and ending pos.
	my $seq = $_[0];
	
	#STARTING POS. OF CDR3:
	my $left_area = substr($seq, 92, 11); #look from pos. 92 - 102 of seq (0-based)
	my $la_i; my $left_cdr = -1;
	if($left_area =~ /(Y.C)/)
	{#if we found 'YXC', find its index
		$la_i = index($left_area, $1);
		$la_i += 2; #'C' is our marker, so advance past 'YX'
	}
	else { $la_i = index($left_area, 'C'); } #didn't find 'YXC', look for 'C'
	if($la_i >= 0) { $left_cdr = $la_i + 92 + 3; } #CDR3 starts at 'C' + 3 (add 92 to put it back in the full sequence)
	
	#ENDING POS. OF CDR3:
	my $n = length($seq)-1; my $n1 = $n-14;
	my $subtract_amount = 1; 
	my $right_area = substr($seq, $n1, 11); #look from pos. n-14 - n-4 of seq (n = last index of seq)
	my $ra_i; my $right_cdr = -1;
	if(($ra_i = index($right_area, 'WGQ')) < 0)
	{#didn't find 'WGQ', look for 'WG'
		if(($ra_i = index($right_area, 'WG')) < 0)
		{#didn't find 'WG', look for 'W'
			if(($ra_i = index($right_area, 'W')) < 0)
			{#didn't find 'W', look for 'GQ'
				$ra_i = index($right_area, 'GQ');
				if($ra_i >= 0) { $subtract_amount = 2; } #if 'GQ' found, CDR3 ends at 'G' - 2 
			}
		}
	} 
	if($ra_i >= 0) { $right_cdr = $ra_i + $n1 - $subtract_amount; } #CDR3 ends at 'W' - 1 (or 'G' - 2) (add n-14 to put it back in the full sequence)
	
	
	#check if st/end found and if not follow rules: 
	if($left_cdr == -1 && $right_cdr == -1) { $left_cdr = $n-21; $right_cdr = $n-10; }
	elsif($left_cdr == -1) { $left_cdr = $right_cdr - 11; }
	elsif($right_cdr == -1) { $right_cdr = ($left_cdr + 11) <= $n ? ($left_cdr + 11) : $n; }
	
	if($left_cdr > $right_cdr)
	{
		#print "Alert! left_cdr ($left_cdr) > right_cdr ($right_cdr) in find_cdr3: $seq\n";
		$left_cdr = $n-1;
		$right_cdr = $n;
	} 
	
	$_[1] = $left_cdr;
	$_[2] = $right_cdr;
	
	return 1;
}

sub compare_cdrs
{#return true if the 2 sequences differ by 0 or 1 AA
	my @cdr1 = split '', $_[0];
	my @cdr2 = split '', $_[1];
	
	if(abs($#cdr2-$#cdr1) > 1)  { return 0; }
	
	my $num_diff = $#cdr2 > $#cdr1 ? ($#cdr2-$#cdr1) : 0;
	
	my $i = 0;
	foreach (@cdr1)
	{
		if(($i > $#cdr2) || ($cdr2[$i] ne $_)) { $num_diff++; }
		if($num_diff > 1) { return 0; }
		$i++;
	}
	
	#if($num_diff > 1) { return 0; }
	#else { return 1; }
	return 1;
}

sub consolidate_groups
{#adds group $j to group $i, if a matching sequence in both groups is found ($nanobody_group[$i] should be defined)
 #then consolidates all previous groups (back to $i) (in case adding group $j caused additional matches w/ previous groups)
	my $i = shift;
	my $j = shift;
	
	#if($nanobody_group[$i])
	if($j <= $i) { return 0; }
	
	my $consolidate = 0;
	if(!defined ${$nanobody_group[$j]}{-1} && !defined ${$nanobody_group[$i]}{-1}) #make sure these groups weren't removed already
	{
		if(!defined ${$nanobody_group[$i]}{$j})
		{
			foreach (keys %{$nanobody_group[$j]})
			{
				if(defined ${$nanobody_group[$i]}{$_})
				{ $consolidate = 1; last; }
			}
		}
		else { $consolidate = 1; }
	}
	
	if($consolidate) 
	{#put $j and all in its group into the $i group
	 #remove $j group
		${$nanobody_group[$i]}{$j} = 1;
		foreach (keys %{$nanobody_group[$j]})
		{
			${$nanobody_group[$i]}{$_} = 1;
		}
		%{$nanobody_group[$j]} = ();
		${$nanobody_group[$j]}{-1} = 1; #mark this group as removed (different than a sequence that matches no other sequences)
		
		
	}
	consolidate_groups($i, $j-1);
}

sub close_OUT_ALL
{
	my $group_div = shift;
	my $next_num = shift;
	
	if($group_div) { print OUT_ALL qq!</div>!; }
	print OUT_ALL "</pre>";
	
	if ($next_num > 1)
	{
		my $cur_num = $next_num-1;
		#put link to go to next page
		print OUT_ALL qq(<br><br>page $cur_num | <a href="$fileroot.$next_num.cdr_coverage.html">next page</a>);
	}
	
	close(OUT_ALL);
}

sub open_OUT_ALL
{
	my $num = shift;
	
	#if($num == 1) { open(OUT_ALL,">$filedir/$fileroot.cdr_coverage.html") || die "Could not open for writing: $filedir/$fileroot.cdr_coverage.html\n"; }
	#else { open(OUT_ALL,">$filedir/$fileroot.$num.cdr_coverage.html") || die "Could not open for writing: $filedir/$fileroot.$num.cdr_coverage.html\n"; }
	if(!open(OUT_ALL,">$filedir/$fileroot.$num.cdr_coverage.html"))
	{
		print LOG "Could not open for writing: $filedir/$fileroot.$num.cdr_coverage.html ($!)\n";
		close(LOG);
		exit(1);
	}
	
	print OUT_ALL <<JSCRIPT;
	<SCRIPT LANGUAGE="JavaScript">
	// Row Hide function.
	function ec(tId, clickIcon)
	{
		dstyle = document.getElementById(tId).style.display;
		if (dstyle == "none")
		{
			document.getElementById(tId).style.display = "";
			document.getElementById(clickIcon).src = "$img_source_minus";
			document.getElementById(clickIcon).alt = "-";
		}
		else
		{
			document.getElementById(tId).style.display = "none";
			document.getElementById(clickIcon).src = "$img_source_plus";
			document.getElementById(clickIcon).alt = "+";
		}
	}
	</SCRIPT>
JSCRIPT
	print OUT_ALL "<pre>\n";
}

