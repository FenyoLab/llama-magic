#!/usr/local/bin/perl
#
#    parse_xtandem_llama.pl
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

my $error=0;
my $dir="";
my $expect_threshold="";
if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="C:\\NCDIR\\Llama\\results\\28\\32\\tandem\\results\\"; }
if ($ARGV[1]=~/\w/) { $expect_threshold=$ARGV[1];} else { $expect_threshold=0.1; } #1e-2; }

sub log10
{
	my $n = shift;
	return log($n)/log(10);
}

if ($error==0)
{
	$dir=~s/\\/\//g;
	#open (LOG,qq!>$dir.$expect_threshold.log!) || die "Could not open output\n";
	my $line="";
	if (opendir(DIR,"$dir"))
	{
		my @allfiles=readdir DIR;
		closedir DIR;
		foreach my $filename (@allfiles)
		{
			if ($filename=~/\.xml$/i)
			{
				open (IN,"$dir/$filename") || die "Could not open $dir/$filename\n";
				open (OUT,qq!>$dir/$filename.peptide_list.$expect_threshold.txt!) || die "Could not open output\n";
				#print OUT qq!spectrum\texpect\tlog(e)\tsequence\tprotein\n!;
				print OUT qq!sequence\tlog(e)\tprotein_uid\tdomain_id\tspectrum\n!;
				#my $peptide_proteins="";
				
				#my $reversed=1;
				#my $pep="";
				#my $expect="";
				#my $id = "";
				#my $uid = "";
				
				my @reversed_list = ();
				my @pep_list = ();
				my @expect_list = ();
				my @id_list = ();
				my @uid_list = ();
				
				while ($line=<IN>)
				{ #<protein expect="-282.5" id="368.1" uid="10753" label="sp|CAH2_HUMAN|" sumI="9.91" >
					#if ($line =~ /NTLYLQMNSMK/) #/"SGTWWYQR/)
					#{
					#	print $line;
					#}
					
					if ($line=~/^\<protein\s+.*uid="([^\"]+)" label="([^\"]+)"/)
					{
						#$uid = $1;
						push(@uid_list, $1);
						
						my $protein_name=$2;
						#my $protein=$protein_name;
						#$protein=~s/^(\S+)\s.*$/$2/;
						if ($protein_name!~/\:reversed$/)
						{
							#$reversed=0;
							push(@reversed_list, 0);
						}
						else { push(@reversed_list, 1); }
						#$peptide_proteins.="#$protein#";
					}
					if ($line=~/^\<domain\s+id="([0-9\.edED\+\-]+)".*expect="([0-9\.edED\+\-]+)".*mh="([0-9\.edED\+\-]+)".*delta="([0-9\.edED\+\-]+)".*seq="([A-Z]+)"/)
					{
						if ($#id_list < $#uid_list)
						{#skip additional domains for the protein (only take first)
							#$id = $1;
							push(@id_list, $1);
							
							my $expect_=$2;
							my $pep_=$5;
							#$pep_=~tr/L/I/;
							#$pep=~tr/L/I/;
							#if ($expect !~ /\w/ or $expect_ <= $expect)
							#{
							#	$expect=$expect_;
							#	$pep=$pep_;
							#}
							push(@expect_list, $expect_);
							push(@pep_list, $pep_);
						}
					}
					if($line=~/<note label=\"Description\">(.+?)<\/note>/)	
					{
						my $title=$1;
						if ($#uid_list != $#reversed_list or $#uid_list != $#expect_list or
						    $#uid_list != $#id_list or $#uid_list != $#pep_list)
						{
							print "Alert!  Error in reading/parsing the xtandem input file!\n";
							
						}
						
						for(my $i=0; $i<=$#uid_list; $i++)
						{
							if ($reversed_list[$i]==0 and $expect_list[$i]<$expect_threshold)
							{
								my $log_e = sprintf("%.1f", log10($expect_list[$i]));
								#print OUT qq!$title\t$expect\t$log_e\t$pep\t$peptide_proteins\n!;
								print OUT qq!$pep_list[$i]\t$log_e\t$uid_list[$i]\t$id_list[$i]\t$title\n!;
							}
						}
						undef @reversed_list; @reversed_list = ();
						undef @pep_list; @pep_list = ();
						undef @expect_list; @expect_list = ();
						undef @id_list; @id_list = ();
						undef @uid_list; @uid_list = ();
						
						#$reversed=1;
						#$pep="";
						#$expect="";
						#$id = "";
						#$uid = "";
						#$peptide_proteins="";
					}
				}	
				close(IN);
				close(OUT);	
			}
		}
	}
	#close(LOG);
}