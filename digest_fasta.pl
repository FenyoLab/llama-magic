#!/usr/local/bin/perl 

#    digest_fasta.pl
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
#require "./masses_and_fragments.pl";

# reads each file in the directory - either the first argument on the command line, or the current directory
#for each fasta file in the directory, the name/sequence is read in and then the protein is digested with 
#trypsin to get the resulting peptides - 
#a single output file is created - "all_predigested.fasta" that contains, 
#for the description line: the peptide sequence, followed by a number representing the number of proteins that 
#  resulted in that peptide when digested with trypsin
#for the sequence line: the peptide sequence

use warnings;
use strict;

my $USE_TAIL = 0;

my $dir="";
my $incompletes=0;

if ($ARGV[0] && $ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="C:/temp/temp/temp"; } #{ $dir="C:/NCDIR/Llama/results/36/protein"; } #"C:\\NCDIR\\Llama\\FASTQ\\Merged\\10_07_2013\\4761_merged_iden0.9_notrim"; } 
if ($ARGV[1] && $ARGV[1]=~/\w/) { $incompletes=$ARGV[1];} else { $incompletes=1; }

my $total_pep_count=0;
my %PEP=();
my %PEP_proteins=();
my %PEP_proteins_count=();
my $files_count=0;
my $peptides_count=0;
my $proteins_count=0;
my %proteins_count=();
my $line="";

if (!opendir(DIR,"$dir")) { print "Error reading $dir\n"; exit(1); }

my @allfiles=readdir DIR;
closedir DIR;
my $count_fasta = 0;
foreach my $filename (@allfiles)
{
	if ($filename!~/\.fas?t?a?$/i) { next; } #can be *.fa, *.fas, *.fast, *.fasta, case insensitive
	
	if (!open (IN,"$dir/$filename")) { print "Error opening $dir/$filename.\n"; next; }
	$count_fasta++;
	#print qq!$filename\n!;
	my $name="";
	my $sequence="";
	my $protein_count = 0;
	while ($line=<IN>)
	{
		chomp($line);
		if ($line =~ s/^>//)
		{#if the current line starts with a '>', then it is the description line, remove the '>' and input the 'name'
			
			my $name_=$line;
			if ($name_ =~ s/^\"//)
			{#if the line begins with quotes (")
			 #remove quotes at beginning and end
			 #replace any non-word character (not letters or numbers) with '_'
				$name_ =~ s/\"$//;
				$name_ =~ s/[^\w]/_/g; 
			}
			else
			{#line does not begin with quotes ("), remove anything after (and including) the first whitespace character 
				$name_ =~ s/\s.*$//;
			}
			if ($name =~ /\w/ and $sequence =~ /\w/)
			{#both name and sequence were inputted, add the protein to the count
				$proteins_count++;
				$proteins_count{$filename}++;
				#print "$proteins_count. $name\n";
				
				#add tail sequence to protein so that we don't miss c terminus peptides:
				if($USE_TAIL) { $sequence .= "SEPKIPQPQPKPQ"; }
				
				#digest the protein with trypsin, the function fills the %PEP hash with the resulting peptides
				%PEP=();
				DigestTrypsin($name,$sequence,$incompletes);
				
				if($protein_count % 100000 == 0) { print "$protein_count\n"; } #if($protein_count >= 10000) { last; } }
				$protein_count++;
				
				foreach my $peptide (keys %PEP)
				{#we want to keep a count of the number of proteins that contained a particular peptide when digested 
				 #with trypsin - a list of each 'name' of the proteins is stored in the hash %PEP_proteins and the number 
				 #of proteins for a given peptide is stored in the hash %PEP_proteins_count
					if (length($peptide)>6)
					{#disregard peptides that are too short
						#if (!$PEP_proteins{$peptide} || ($PEP_proteins{$peptide} !~ /#$name#/))
						#{
						#	$PEP_proteins_count{$peptide}++;
						#	$PEP_proteins{$peptide} .= qq!#$name#!;
						#}
						#else
						#{
						#	print "$name $peptide\n";
						#}
						$PEP_proteins_count{$peptide}++;
					}
				}
				
			}
			$name=$name_;
			$sequence="";
		}
		else
		{#it is a line of the sequence
			$sequence.="$line";
		}
		
		
	
	}	
	
	#the last protein left in the file, do the same as above
	if ($name =~ /\w/ and $sequence =~ /\w/)
	{
		$proteins_count++;
		$proteins_count{$filename}++;
		#print "$proteins_count. $name\n";
		
		#add tail sequence to protein so that we don't miss c terminus peptides:
		if($USE_TAIL) { $sequence .= "SEPKIPQPQPKPQ"; }
				
		%PEP=();
		DigestTrypsin($name,$sequence,$incompletes);
		foreach my $peptide (keys %PEP)
		{
			#print qq!$peptide\n!;
			if (length($peptide)>6)
			{
				#if (!$PEP_proteins{$peptide} || ($PEP_proteins{$peptide} !~ /#$name#/))
				#{
				#	$PEP_proteins_count{$peptide}++;
				#	$PEP_proteins{$peptide} .= qq!#$name#!;
				#}
				$PEP_proteins_count{$peptide}++;
			}
		}
	}
	close(IN);
	$files_count++;
	print qq!$files_count. $filename $proteins_count{$filename}\n!;
	
	
}

if ($count_fasta == 0)
{
	print "Warning: No fasta files processed!\n";
	exit(1);
}

if (open (OUT,">$dir/all_predigested.fasta"))
{
	foreach my $peptide (keys %PEP_proteins_count) #(keys %PEP_proteins)
	{
		print OUT qq!>$peptide $PEP_proteins_count{$peptide}\n$peptide\n!;
	}
	close(OUT);
}
else
{
	print "Error creating $dir/all_predigested.fasta\n";
	exit(1);
}

exit(0);

sub DigestTrypsin
{#fills the %PEP hash with the peptides resulting from digesting the sequence with trypsin
 #arg1 - name, arg2 = sequence, arg3 = # of incompletes
	my $name = shift();
	my $seq = shift();
	my $incompletes = shift();

	my $temp=$seq;
	my @pep=();
	my @start=();
	my @end=();
	my $aa="";
	my $aa_="";
	my $i=0;

	for($i=0;$i<=$incompletes;$i++)
	{
		$start[$i]=0;
		$end[$i]=-1;
		#$pep[$i]="[";
	}
	my $aa_count=0;
	while ($temp =~ s/^\s*([A-Z\*])//)
	{
		$aa="\U$1";
		$aa=~s/I/L/g;
		if ( (($aa_=~/R/ or $aa_=~/K/) and $aa!~/P/) or $aa_=~/\*/)
		{
			for($i=0;$i<=$incompletes;$i++)
			{
				$PEP{"$pep[$i]"}=1;
				$pep[$i]=$pep[$i+1];
				$start[$i]=$start[$i+1];
				$end[$i]=$end[$i+1];
			}
			$start[$incompletes]=$aa_count;
			$end[$incompletes]=$aa_count-1;
		}
		for($i=0;$i<=$incompletes;$i++)
		{
			if ($aa!~/\*/) { $pep[$i].=$aa; }
			$end[$i]++;
		}
		$aa_=$aa;
		$aa_count++;
	}
	for($i=0;$i<=$incompletes;$i++)
	{
		$PEP{"$pep[$i]"}=1;
	}
}


