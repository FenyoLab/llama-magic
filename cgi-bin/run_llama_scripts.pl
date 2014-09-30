#!C:/Perl/bin/perl 

#    run_llama_scripts.pl
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

#runs the scripts sent in by command line arguments for (CGI web script) llama_magic.pl
#options are:
#db_scripts (runs both dnatoprot_longest_nr.pl and digest_fasta.pl)
#search_and_map_scripts (runs xtandem (command line/xml interface), parse_xtandem.pl, filter_tandem_results.pl, and map_peptides_to_proteins.pl)

#use File::Copy;


use warnings;
use strict;

my $expect_threshold = '0.1';

my $function = $ARGV[0];

if ($function eq 'db_scripts')
{
    my $results_dir = $ARGV[1];
    
    #create status txt file
    open(OUT, ">$results_dir/protein/status.txt") or die "Failed to create status file: $results_dir/protein/status.txt";
    
    #run dnatoprot
    print OUT qq!Calling: "../dnatoprot_longest_nr.pl" "$results_dir/dna" "$results_dir/protein"\n!;
    my $cmd_out = `"../dnatoprot_longest_nr.pl" "$results_dir/dna" "$results_dir/protein" 2>&1`;
    if ( $? == -1 )
    {
        print OUT "ERROR: command failed (dnatoprot_longest_nr.pl): $!\n";
    }
    elsif($? >> 8 != 0) # exit value not 0, indicates error...
    {
        printf OUT "ERROR: command (dnatoprot_longest_nr.pl) exited with value %d\n", $? >> 8;
        print OUT "$cmd_out\n";
    }
    else
    {
        print OUT "Success: $cmd_out\n";
        
        #Success! add vh_prot.fas sequences to the longest_nr file:
        #CORRECTION - this file not needed
        #if(!open(PROTEINS_IN, "../vh_prot.FAS"))
        #{
        #    print OUT "ERROR: Opening of '../vh_prot.FAS' for reading failed: $!\n";
        #}
        #else
        #{
        #    if(!open(PROTEINS_OUT, ">>$results_dir/protein/longest_nr.fasta"))
        #    {
        #        print OUT "ERROR: Opening of '$results_dir/protein/longest_nr.fasta' for appending failed: $!\n";
        #    }
        #    else
        #    {
        #        print PROTEINS_OUT "\n"; #just in case
        #        while(<PROTEINS_IN>)
        #        {
        #            print PROTEINS_OUT $_;
        #        }
        #        
        #        close(PROTEINS_IN);
        #        close(PROTEINS_OUT);
        #        print OUT "Success: Added vh_prot.FAS to longest_nr.fasta.\n";
            
                #success! run digest_fasta
                print OUT qq!Calling: "../digest_fasta.pl" "$results_dir/protein" "1"\n!;
                $cmd_out = `"../digest_fasta.pl" "$results_dir/protein" "1" 2>&1`;
                if ( $? == -1 )
                {
                    print OUT "ERROR: command failed (digest_fasta.pl): $!\n";
                }
                elsif($? >> 8 != 0) # exit value not 0, indicates error...
                {
                    printf OUT "ERROR: command (digest_fasta.pl) exited with value %d", $? >> 8;
                    print OUT "$cmd_out\n";
                }
                else
                {
                    print OUT "Success: $cmd_out\n";
                    
                    #success!  rename all_predigested.fasta to longest_nr_predigested.fasta
                    if (!rename "$results_dir/protein/all_predigested.fasta", "$results_dir/protein/longest_nr_predigested.fasta")
                    {
                        print OUT "ERROR: Could not rename '$results_dir/protein/all_predigested.fasta' to '$results_dir/protein/longest_nr_predigested.fasta'.\n";
                    }
                }
            #}
        #}
        
    }
    #write to status file and close
    print OUT "DONE\n";
    close(OUT); 
}
elsif($function eq 'search_and_map_scripts')
{
    my $results_dir = $ARGV[1];
    my $taxon_dir = $ARGV[2];
    my $parent_err = $ARGV[3];
    my $fragment_err = $ARGV[4];
    my $cgi_results_dir = $ARGV[5];
    my $show_score = $ARGV[6];
    
    #create status txt file
    open(OUT, ">$results_dir/status.txt") or die "Failed to create status file: $results_dir/status.txt";
    
    #create xtandem input.xml
    if(open(XML_OUT, ">$results_dir/tandem/input.xml"))
    {
        print XML_OUT <<XMLTEXT;
<?xml version="1.0"?>
<bioml>
	<note type="input" label="spectrum, fragment monoisotopic mass error">$fragment_err</note>
	<note type="input" label="spectrum, parent monoisotopic mass error plus">$parent_err</note>
	<note type="input" label="spectrum, parent monoisotopic mass error minus">$parent_err</note>
	<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>
	<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>

	<note type="input" label="list path, default parameters">../default_input.xml</note>
	<note type="input" label="list path, taxonomy information">$taxon_dir/taxonomy.xml</note>

	<note type="input" label="protein, taxon">llama</note>
XMLTEXT

        #get list of mgf files and add a line in input.xml for each one
        my @file_names = <$results_dir/mgf/*.mgf>; 
        foreach my $cur_file (@file_names)
	{
            print XML_OUT qq(\n\t<note type="input" label="spectrum, path">$cur_file</note>);
        }
        
        print XML_OUT <<XMLTEXT;
	<note type="input" label="output, path">$results_dir/tandem/results/output.xml</note>
	<note type="input" label="output, results">valid</note>
</bioml>
XMLTEXT

        close(XML_OUT);
    
        #run xtandem with uploaded mgf files
        print OUT qq!Calling: "../xtandem.exe" "$results_dir/tandem/input.xml"\n!;
        my $cmd_out = `"../tandem.exe" "$results_dir/tandem/input.xml" 2>&1`;
        if ( $? == -1 )
        {
            print OUT "ERROR: command failed (tandem.exe): $!\n";
        }
        elsif($? >> 8 != 0) # exit value not 0, indicates error...
        {
            printf OUT "ERROR: command (tandem.exe) exited with value %d\n", $? >> 8;
            print OUT "$cmd_out\n";
        }
        else
        {
            #parse xtandem to tab separated txt file, extract relevant columns
            print OUT qq!Calling: "../parse_xtandem_llama.pl" "$results_dir/tandem/results" "$expect_threshold"\n!;
            my $cmd_out = `"../parse_xtandem_llama.pl" "$results_dir/tandem/results" "$expect_threshold" 2>&1`;
            if ( $? == -1 )
            {
                print OUT "ERROR: command failed (parse_xtandem_llama.pl): $!\n";
            }
            elsif($? >> 8 != 0) # exit value not 0, indicates error...
            {
                printf OUT "ERROR: command (parse_xtandem_llama.pl) exited with value %d\n", $? >> 8;
                print OUT "$cmd_out\n";
            }
            else
            {
                #get name of txt file created above:
                #my @files = <$results_dir/tandem/results/*.txt>; 
                #my $pep_file = $files[0]; #should be exactly 1 txt file in directory
                
                my $pep_file = "output.xml.peptide_list.$expect_threshold.txt";
                
                #run map_peptides_to_proteins - output is candidate list html file
                
                print OUT qq!Calling: "../map_peptides_to_proteins.pl" "$results_dir/tandem/results/$pep_file" "$taxon_dir/protein/longest_nr.fasta" "$cgi_results_dir/tandem/results/output.xml" "$show_score"\n!;
                my $cmd_out = `"../map_peptides_to_proteins.pl" "$results_dir/tandem/results/$pep_file" "$taxon_dir/protein/longest_nr.fasta" "$cgi_results_dir/tandem/results/output.xml" "$show_score" 2>&1`;
                if ( $? == -1 )
                {
                    print OUT "ERROR: command failed (map_peptides_to_proteins.pl): $!\n";
                }
                elsif($? >> 8 != 0) # exit value not 0, indicates error...
                {
                    printf OUT "ERROR: command (map_peptides_to_proteins.pl) exited with value %d\n", $? >> 8;
                    print OUT "$cmd_out\n";
                }   
            }
            
        }
    }
    else { print OUT "ERROR: Failed to create xtandem input file: $results_dir/tandem/input.xml"; }     
    
    #write to status file and close
    print OUT "DONE\n";
    close(OUT); 
}


