#!C:/Perl/bin/perl -d
##
##
## peptide.pl
## Copyright (C) 2003-2004 Ronald C Beavis, all rights reserved
## The Global Proteome Machine 
## This software is a component of the X! proteomics software
## development project
##
## Use of this software governed by the Artistic license,
## as reproduced at http://www.opensource.org/licenses/artistic-license.php
##
## Version 2004.01.07
## Version 2004.02.01
## Version 2004.02.02
## Version 2004.03.01
## Version 2004.08.01
## Version 2004.08.01a - new style
## Version 2004.09.20
## Version 2004.09.27 - solves issue in draw_fragment() when the delta masses are set. arises from
##						situations where two ions have a mass difference of m_fAmmonia. This causes
##						the delta loop to be exited prematurely. now the loop continues through each
##						possible ion mass. Also, if a b and y ion have very similar masses, both are shown in 'fragmentogram'.
## Version 2004.09.30 - get n and c terminal mass changes from xml. Hover states for delta diagram added to show mass of ion.
## Version 2004.10.04 - add check for expect > 0 before doing log()
## Version 2004.10.15 - add call to get_root()
## Version 2004.10.29 - add proex var; the minimum protein expect. changed print statements to print qq(
## Version 2004.11.09 - adds in the mz and intensity values to the validate link to facilitate angle filtering
## Version 2004.11.15 - check for .gz result, extract and process if exists, otherwise use original
## Version 2004.12.03 - added link to peptide atlas
## Version 2005.01.01 - altered display to account for neutral losses
## Version 2005.02.15 - changed mass accuracy & values to agree with X! Tandem v. 2005.02.01
## Version 2005.02.22 - 1. get residue, H2O and NH3 mass values from output.xml if 
##						<group label="residue mass parameters" type="parameters"> exists in output.xml file.
##						this is done in set_aa() function, which has been moved to common.pl. It takes
##						the output.xml file as a parameter and returns a hash ref containing the mass values.
##						Also, use $m_fWater as value for $m_fY (value lost on c-terminal cleavage)
##						2. If delta attribute does not exist, try to get other domain info. If this still fails,
##						report an error.
## Version 2005.02.25 - changes to is_special() to check for phosphorylation. If is_special() also add mass and intensity 
##						to marked mass array
##
## Version 2005.03.10 - 1. changes to fragment-o-gram, histogram and delta display to show mz for +2 and +3 ions
##						2. add hover states for show/hide of matched/unmatched ions on histogram
##
## Version 2005.03.14 - 1. add -H20 to all levels of display using aqua (b ion) and purple (y ion) for colors.
##						2. removed space from beginning and end of mass and intensity lines.
##						3. fix alignment on z column
##
## Version 2005.03.28 - 1. update parameters for SystemsBiology link
## Version 2005.08.26 - 1. added detection of +2 y ions with P terminus
## Version 2005.09.19 - 1. added detection of methyl-sulphone neutral loss
## Version 2005.10.19 - 1. added detection of +2 y ions with any terminus
## Version 2006.04.04 - fixed problem with version 3.03 for SVG ActiveX control
##
## Version 2006.06.14 - javascript fix for svg requiring click to activate
## Version 2006.07.18 - addition of Pepseeker database link
## Version 2007.02.28 - corrected links to Pepseeker and PeptideAtlas
## Version 2007.06.13 - added SwedCAD link
## Version 2007.08.02 - detection of neutral loss of urea adduct
## Version 2007.08.07 - removed "prompt" fragmentation from S/T residues unless specified
## Version 2007.12.06 - added AJAXed peptide count for reference
## Version 2008.01.27 - correct errors displaying unusual mutated sequences
## Version 2008.04.10 - added MRM link
## Version 2008.07.15 - added mass range value to spectrum validation form
## Version 2008.07.15 - corrected error in MRM link
## Version 2008.07.31 - altered MRM link
## version 2008.10.09 - added AJAX-fetched omega score
## version 2009.05.21 - modified omega display control rendering to check if the label appears 
##	in the mrm database, instead of whether or not the label matches a list of regexes.

## peptide.pl creates an svg image that represents the peptide by using
## the peaks assigned by tandem. In the same svg image it creates an error
## graph and a pictorial representation of the sequence showing the assigned b and y ions.
## At the end of the page it displays all theoretical and assigned ions for the peptide
## in a color coded table layout.

##  parameters (CGI): 
##		label - the "accession number" eg: ENSP00000244534 or gi|123456
##		homolog - the protein uid to search for
##		path - the path to the current result file
##		id - the id of the <domain>
##		uid - the uid used in the "get annotation" link 
##	called by: protein.pl
##	required by: none
##	calls/links to: none


use strict;
use CGI qw(:all);
use CGI::Carp qw(fatalsToBrowser);
use LWP::UserAgent;

require "./common.pl";
require "./defines.pl";
my $g_server = get_server_name();

my $file_version = "peptide.pl, v. 2011.12.21";


my $gvalue = GetGsite('gpmdb_url');
my $wvalue = GetGsite('wiki_url');
my $mvalue= GetGsite('mrm_url');

my $cgi = CGI->new();

my $use_png = 0;
##set to one if domain info is incomplete
my $error = 0;

my $path = $cgi->param('path');
if ($path!~/^\/\//) { $path = get_root() . $cgi->param('path'); }
my $url = $cgi->param('path');
my $label = $cgi->param('label');
my $uid_input = $cgi->param('uid');
my $homolog = $cgi->param('homolog');
my $id = $cgi->param('id');
my $proex = $cgi->param('proex');
my ($gpm) = $url =~ /(GPM[0-9]+)/;
my $ltype = $cgi->param('ltype');

my $bOmega = 0;
my $ua = LWP::UserAgent->new;
my $usu = GetGsite('gpmdb_user');
my $usp = GetGsite('gpmdb_pwd');

my $urla = 'http://';
if($usu and $usp)	{
	$urla .= "$usu:$usp@";
}

$urla .= $g_server . "/thegpm-cgi/request_server_mrm.pl?target=mrm_pro_check&data=$label";
#my $resp = $ua->get($urla); #LLAMA CHANGES

my $pro_check_text = "";
#my $pro_check_text = $resp->content; #LLAMA CHANGES
#$pro_check_text =~ s/<[^>]+>//g; #LLAMA CHANGES

my @spec_matched;
my @spec_unmatched;
my @spec_total;
my %spec_stats;
my $spec_gpm = $gpm;
my $spec_id = $id;

my $mgf = $cgi->param('unmatched');	
if($mgf)	{
	generate_mgf();
	exit;
}

if ($pro_check_text =~ /\=[1-9]/) { # at least one record exists for this protein

	$bOmega = 1;
}

PrintHeader("peptide model: $id of $label");

## set arrays for sequence calculations
my $m_pfAaMass=set_aa($path);
my $m_fWater = $m_pfAaMass->{'H2O'};
my $m_fAmmonia = $m_pfAaMass->{'NH3'};
my @m_pSequence;
my @m_pfMods;
my $m_lSequence;
my @m_pfY;
my @m_pfI;
my @m_strI;
my @m_pfB;
my @m_pfYd;
my @m_pfBd;
my $m_sigma;

##+2 and +3 charged deltas
my @m_pfYdx2;
my @m_pfBdx2;
my @m_pfYdx3;
my @m_pfBdx3;
##

my @m_pfZ;
my @m_pfZ1;
my @m_pfC;
my @m_pfZd;
my @m_pfZ1d;
my @m_pfCd;
my @m_pfZdx2;
my @m_pfZ1dx2;
my @m_pfCdx2;
my @m_pfZdx3;
my @m_pfZ1dx3;
my @m_pfCdx3;
my @itable;

my $g_cid = 0;
my $g_etd = 0;

## set defaults for sequence calculations
my $m_fProton = 1.007276;
my $m_fHydrogen = 1.007825035;

my $m_fN = 0.0;
my $m_fC = 17.026549105;
my $m_fZ = $m_fWater - $m_fAmmonia + $m_fProton;
my $m_fA = -27.99491463;
my $m_fB = 0.0;
my $m_fY = $m_fWater;
my $m_fX = 43.990384;
my $m_fCleaveC = 17.002739665;
my $m_fCleaveN = 1.007825035;
my $m_fCleaveCdefault = 17.002739665;
my $m_fCleaveNdefault = 1.007825035;
my $m_fError = 0.5;
my $m_pErrorType = "Daltons";
my $m_fMassMax = 0.0;

my $prompt_phospho = 0;

my $line;
my $e;
my $feature;
my $domain;
my $start;
my $end;
my $delta;
my $mass;
my $expect;
my $bI = 1;
my $intensity=0;
my $activation = 0;
my $charge;
my @res;
my @res_pos;
my @res_mod;
my @res_mut;
my @res_prompt;
my $mass_line;
my $int_line;
my $note_line = "";
my $id_v = $id;
my $annotation;
my $ann;
$id_v =~ s/\./\\\./g;

$/ = "<\/group><\/group>"; # we can reset the record separator from newline (default) to whatever we want.  This brings in
                         # one model at a time as a single string.
my $e = $id;
$e =~  s/\..+/\.spectrum/;

my %m_IonTypes;
my $s_pre;
my $s_post;
if(not open(INPUT,"<$path")){
	if(open(INPUT,"<$path.gz"))	{
		close(INPUT);
		system("gzip -d $path.gz");
		open(INPUT,"<$path");
	}
	else	{
		PrintHeader("");
		print qq(
			<table><tr><td>
			<a href="/index.html"><img src="/llama-magic-html/gpm.png" border="0"></a>
			</td><td>&nbsp;&nbsp;</td><td valign="middle" width="400">
			<i>$label</i><BR><BR>
			The model file <b>$path</b> <BR>
			could not be found in the archive.</td></tr>
			</table></div></div></body></html>\n);
		exit;
	}
}
my $gpm_show = 0;
my %gpm_value;
my %accessions;
my @fixed_mods;
my @variable_mods;
while(<INPUT>)	{
	if(/\<group/ and /label\=\"input parameters\"/)	{
	  	if (/spectrum, fragment monoisotopic mass error units\">(\w+)</){ 
  			$m_pErrorType = $1;
 	 	}
 	 	if (/spectrum, fragment monoisotopic mass error\">(.+?)</){
  			$m_fError = $1;
 	 	}
  		## rc - get n and c terminal mods
 	 	if (/protein, cleavage N-terminal mass change\">(.+?)</){
 			$m_fCleaveN = $1;
 	 	}
 		if (/protein, cleavage C-terminal mass change\">(.+?)</){
 			$m_fCleaveC = $1;
 		}
		if(/label\=\"scoring, a ions\"/)	{
			if(/yes/)	{
				$m_IonTypes{"a"} = 1;
			}
		}
		if(/label\=\"scoring, b ions\"/)	{
			if(/yes/)	{
				$m_IonTypes{"b"} = 1;
			}
		}
		if(/label\=\"scoring, c ions\"/)	{
			if(/yes/)	{
				$m_IonTypes{"c"} = 1;
			}
		}
		if(/label\=\"scoring, x ions\"/)	{
			if(/yes/)	{
				$m_IonTypes{"x"} = 1;
			}
		}
		if(/label\=\"scoring, y ions\"/)	{
			if(/yes/)	{
				$m_IonTypes{"y"} = 1;
			}
		}
		if(/label\=\"scoring, z ions\"/)	{
			if(/yes/)	{
				$m_IonTypes{"z"} = 1;
			}
		}
		my @lines = split /\n/,$_;
		my $l;
		my $x;
		foreach $l(@lines)	{
			if($l =~ /\<\/group/)	{
				last;
			}
			if($l =~/\"residue\, modification mass/)	{
				($x) = $l =~ /\>(.+?)\</;
				if(length($x) > 3)	{
					push(@fixed_mods,$x);
				}
			}
			if($l =~/\"residue\, potential modification mass/)	{
				($x) = $l =~ /\>(.+?)\</;
				if(length($x) > 3)	{
					unshift(@variable_mods,$x);
				}
			}
			if($l =~/\"refine\, potential modification mass/)	{
				($x) = $l =~ /\>(.+?)\</;
				if(length($x) > 3)	{
					push(@variable_mods,$x);
				}
			}
		}
		my (@notes) = $_ =~ /\<note type\=\"input\" label\=\"gpmdb, .+?\>([^\<]+)/gsi;
		my (@names) = $_ =~ /\<note type\=\"input\" label\=\"gpmdb, ([^\"\']+)/gsi;
		my $n = 0;
		while($n < scalar(@names))	{
			@notes[$n] =~ s/^\s+//;
			@notes[$n] =~ s/\s+$//;
			if(length(@notes[$n]))	{
				$gpm_value{@names[$n]} = @notes[$n];
				$gpm_show++;
			}
			$n++;
		}
		if(!$gpm_value{"name"})	{
			$gpm_value{"name"} = "anonymous";
		}
   	}  
	next unless (/domain id=\"$id_v\"/);  # won't even look at the record unless it contains our id
	my @lines = split /\n/,$_;
	my $l;
	foreach $l(@lines)	{
		if($l =~ /\<protein/ and $l =~ /$label/)	{
			$annotation = get_feature($l,"annotation");
		}
		if($l =~ /\<protein/)	{
			$accessions{get_feature($l,"label")} = get_feature($l,"uid");
		}	
	}
  	($charge) = /\<group.+?z=\"(\d)/;
	($activation) = /act\=\"(\d+)\"/;
	if($bI and $charge ne "")	{
		$intensity = get_feature($_,"sumI");
		if($intensity eq "_")	{
			$bI = 0;
		}
	}

  	if(!(($start,$end,$expect,$mass,$delta,$domain) = /<domain id=\"$id_v\" start=\"(\d+).+?end=\"(\d+).+?expect=\"(.+?)\".+?mh=\"(.+?)\".+?delta=\"(.+?)\".+?seq=\"(.+?)\"/)){
  		($start,$end,$expect,$mass,$domain) = /<domain id=\"$id_v\" start=\"(\d+).+?end=\"(\d+).+?expect=\"(.+?)\".+?mh=\"(.+?)\".+?seq=\"(.+?)\"/;
		if(!($start)){
			$error=1;
		}
	}
	($s_pre,$s_post) = /<domain.+?pre=\"(.+?)\".+?post=\"(.+?)\"/;
	my $type;
	my $at;
	my $mod;
	my $prompt;
	my $pm;
	my $aa;
	my $saved = $_;

	my $aastring = "";
	($aastring) = /\<domain id=\"$id_v\".+?>(.+?)\<\/domain>/s; 
	if($aastring ne ""){
		$_ = $aastring;
		my @aas = split /\n/sg;
		foreach $aa(@aas){
			$_ = $aa;
			($type) = /\<aa type=\"(\w)\".+?/g;
			($at) = /\<aa.+?at=\"(\d+)\".+?/g;
			($mod) = /\<aa.+?modified=\"(.+?)\".+?/g;
			($prompt) = /\<aa.+?prompt=\"(.+?)\".+?/g;
			($pm) = /\<aa.+?pm=\"([\w-])\".+?/g;
			if($type eq ""){
				next;
			}
			push(@res,$type);
		    push(@res_pos,$at);
		    if($prompt_phospho and ($type eq 'S' or $type eq 'T'))	{
				if(not $prompt and abs($mod - 80.0) < 1.0)	{
					$prompt = -98.0;
				}
			}
		    push(@res_mod,$mod+$prompt);
		    push(@res_prompt,$prompt);
		    push(@res_mut,$pm);
		}
	}

	$_ = $saved;
	
	# grab the note line (capital 'D'), if exists.
	($note_line) = /<note label=\"Description\">(.+?)<\/note>/g;
	
	# grab mass values
	($mass_line) = /<GAML\:Xdata label=\"$e\".+?numvalues=\"\d+\">(.+?)<\/GAML/s;
	$mass_line =~ s/\s+/ /g;
	$mass_line =~ s/^\s+//;
	$mass_line =~ s/\s+$//;
	
	# grab int values
	($int_line)  = /<GAML\:Ydata label=\"$e\".+?numvalues=\"\d+\">(.+?)<\/GAML/s;
	$int_line =~ s/\s+/ /g;
	$int_line =~ s/^\s+//;
	$int_line =~ s/\s+$//;
}
close(INPUT);
#exit(0);
$/ = "\n";
#		T_Y =	0x01,
#		T_B =	0x02,
#		T_X =	0x04,
#		T_A =	0x08,
#		T_C =	0x10,
#		T_Z =	0x20,
#
if($activation)	{
	%m_IonTypes = ();
	if($activation & 0x01)	{
		$m_IonTypes{"y"} = 1;
	}
	if($activation & 0x02)	{
		$m_IonTypes{"b"} = 1;
	}
	if($activation & 0x010)	{
		$m_IonTypes{"c"} = 1;
	}
	if($activation & 0x04)	{
		$m_IonTypes{"x"} = 1;
	}
	if($activation & 0x08)	{
		$m_IonTypes{"a"} = 1;
	}
	if($activation & 0x20)	{
		$m_IonTypes{"z"} = 1;
	}
}
if($intensity eq "_")	{
	$bI = 0;
}
my @mass_values = split / /,$mass_line;
my @int_values = split / /,$int_line;

StartTop($gpm,"peptide model: $id of $label");

if(scalar(%m_IonTypes) == 0)	{
	$m_IonTypes{"y"} = 1;
	$m_IonTypes{"b"} = 1;
}

my $fixed = $label;
$fixed =~ s/\|//g;
#LLAMA CHANGES
#print qq(
#
#	| <a href="/thegpm-cgi/plist.pl?path=$url&amp;proex=$proex&amp;ltype=$ltype" class="small_link">model</a> | 
#	<a href="/thegpm-cgi/protein.pl?path=$url&uid=$uid_input&amp;label=$label&ltype=$ltype&homolog=$homolog&amp;proex=$proex" class="small_link">protein</a> | 
#	<a href="/thegpm-cgi/homolog.pl?path=$url&label=$label&amp;ltype=$ltype&homolog=$homolog&amp;proex=$proex" class="small_link">homologues</a> | 
#	<a href="$url" class="small_link">XML</a> | 
#	<a href="http://$gvalue/thegpm-cgi/dblist_pep.pl?seq=$domain" class="small_link">gpmDB</a> |
#	<a href="http://$wvalue/wiki/$gpm/$fixed/$domain" class="small_link" target="_WIKI" title="GPM wiki entry for the peptide $domain">wiki</a> | 
#);
if($label =~ /ENSP/){
	print qq(
	  <a href="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide?_tab=3&amp;atlas_build_id=335&amp;searchWithinThis=Peptide Sequence&amp;searchForThis=$domain&amp;query=QUERY" title="Peptide Atlas @ ISB" target="_BLANK">
	Peptide Atlas</a> |
	);
}

my $linfo = translate_label($label,$ltype);
my $info;
if($linfo)	{
	$info = GetInfoNew($label);
}
$info =~ s/\<span.+//si;
$info =~ s/\<br.+//si;
#LLAMA CHANGES
#print qq(<br><br>
#	<table cellspacing='0' cellpadding='0' width='550'><tr class='alt'><td valign='top' align='right' width='150'><b>$linfo</b>:&nbsp;</td><td valign='top' align='left' width='400'>$info</td></tr></table>
#	</td>
#	</tr></table>
#);

print "<script type='text/javascript' src='/llama-magic-html/functions_ajax.js'></script>\n";
print "<script type='text/javascript' src='/llama-magic-html/functions_ajax_mrm.js'></script>\n";
# call ajax functions

my $divctr = 10000; # counter for AJAXed tags
my $divctr_start = $divctr;
my $omegactr = 1000;  # peptide frequency counter
my $omegactr_start = $omegactr;

my $pic_path = $url;
$id =~ s/\..*//;
my $v = int(rand(time|$$)/1000000);
$pic_path ="/results/archive/pics/$id-$v.svg";
my $pic_path_url ="/llama-magic/archive/pics/$id-$v.svg";

print qq(
	<a href="#contributor">Sample information</a><br />
	<table width="726" cellspacing="3" cellpadding="0">
	<tr>
	<td align="center" valign="middle" width="40" title="assignment id number"><b>#</b></td>
	<td align="center" valign="middle" width="50" title="log of E-value for the assignment"><b>log(e)</b></td>
);
if($bI)	{
	print qq(
	<td align="center" valign="middle" width="50" title="log of fragment ion intensities"><b>log(I)</b></td>
	<td align="center" valign="middle" width="60" title="mass of assigned sequence"><b>m+h</b></td>
	<td align="center" valign="middle" width="60" title="mass error for assigned sequence"><b>delta</b></td>
	<td align="center" valign="middle" width="30" title="ratio of measured-to-potential z for assigned sequence"><b>&zeta;</b></td>
	);
}
else	{
	print qq(
	<td align="center" valign="middle" width="80" title="mass of assigned sequence"><b>m+h</b></td>
	<td align="center" valign="middle" width="70" title="mass error for assigned sequence"><b>delta</b></td>
	<td align="center" valign="middle" width="40" title="ratio of measured-to-potential z for assigned sequence"><b>&zeta;</b></td>
	);
}
$v = "http://$g_server$pic_path_url";

if(not($label =~ /\:reversed/))	{
	#LLAMA CHANGES
	print "<td align=\"left\" valign=\"middle\" width=\"400\"><b>sequence</b> "; #| <a href=\"#\" onclick=\"document.forms[1].submit()\" class=\"small_link\">validate</a>";
	#print " | <a href=\"#\" onclick=\"document.forms[2].submit()\" class=\"small_link\">studio</a>";
	print " | <a href=\"#\" onclick=\"document.forms[3].submit()\" class=\"small_link\">mgf</a>";
	my $mz = ($mass - $m_fProton)/$charge + $m_fProton;
	my ($pattern) = $label =~ /(^[A-Z]+)/i;
	if($bOmega)	{
		print " | <a href=\"http://$mvalue/thegpm-cgi/peak_search.pl?pmass=$mz&pmrange=0.1&pattern=$pattern&label=&seq=$domain&unique=&fmrange=0.1&sort=mass&submit=Search\" class=\"small_link\" target=\"_MRM\">mrm</a>";
	}
	#print " | <a href=\"/thegpm-cgi/transform_xml.pl?type=support&id=$id&file=..$url&proex=-1\" title=\"Show additional information supporting this sequence\" class=\"small_link\" target=\"_MRM\">details</a>";
	print " |</td>\n";
}
print("</tr>");
my $a = 0;
my $num = 0;
@m_pSequence = split //,$domain;
$m_lSequence = length($domain);

set_mods();
my $b = 0;
my @seq;
my $aa;
my $value;
my $c;

@seq = split //,$domain;
print("<tr>");
$e = $id;
$e =~ s/\..+//;
print "<td align=\"center\" valign=\"top\">$e</td>";
if($expect > 0){
	$num = log($expect)/2.303;
}
else{
	$num=log(1.0E-10)/2.303;
}
$e = sprintf("%.1f",$num);

print "<td align=\"center\" valign=\"middle\" >$e</td>";

if($bI)	{
	print "<td align=\"center\" valign=\"middle\" >$intensity</td>";
}

my $zeta = zeta($domain,$charge,1);

$s_pre =~ tr/[A-Z]/[a-z]/;
$s_post =~ tr/[A-Z/[a-z]/;

print qq(
	<td align="center" valign="middle" >$mass</td>
	<td align="center" valign="middle" >$delta</td>
	<td align="center" width="30" valign="middle" >$zeta</td>
	<td align="left" valign="middle" >$s_pre<sup>$start</sup>
);

$b = 0;

#LLAMA CHANGES
#print "<a href=\"http://$gvalue/thegpm-cgi/dblist_pep.pl?seq=$domain\" title=\"Observations of this peptide in GPMDB\" target=_BLANK>";
my $js_seq=''; # scalar to hold constructed sequence
for (my $i = 0; $i < scalar@m_pSequence;$i++)	{
	my $print_string;  #initialized but undef
	my $aa = $m_pSequence[$i];
	for(my $c = 0; $c < scalar(@res_pos); $c++)	{
		if($res_pos[$c] == $start + $i)	{
  			if(not($res_mut[$c] eq ""))	{
				$print_string = "<span class=\"mut\">$aa</span>";
  			}
  			else{
				$print_string = "<span class=\"mod\">$aa</span>";
  			}
		}
	}
	print ( (defined$print_string) ? $print_string : $aa);  #fancy "switch" statement.  I could have used if/else
	$js_seq .= $aa; # construct string of residues without markup

}
#print "</a>";
print "<sup>$end</sup>$s_post";
#print "<sup>$end</sup>$s_post&nbsp;&nbsp;<span style='background: #ffffaa;' title='# of observations in GPMDB' id='seqdiv_$divctr' label='$label' charge='$charge' sequence='$js_seq'>(<img src='/llama-magic-html/waiting.gif'>)</span>&nbsp;";

if($bOmega)	{
	print "<span title='peptide frequency in this protein' id='seqdiv_$omegactr' label='$label' charge='$charge' sequence='$js_seq' style='background: #ffffaa;'><img src='/llama-magic-html/waiting.gif'></span>";
	$omegactr++;
}
print "</td>\n";

$divctr++;

undef $e;

my $m_bMeOx = 0;
my $m_bPhospho = 0;
my $m_bUrea = 0;
my $v = 0;
for (my $i = 0; $i < scalar@res; $i++){
	if(defined $e)	{
		$e .= ", ";
	}
	if($res_mut[$i] ne "")	{
		if($res_mod[$i] < 0)	{
			if($res_mut[$i] eq "-")	{
  				$e .= ":p.$res[$i]$res_pos[$i]del ($res_mod[$i])";
			}
			else	{
  				$e .= ":p.$res[$i]$res_pos[$i]$res_mut[$i] ($res_mod[$i])";
			}
		}
		else{
  			$e .= ":p.$res[$i]$res_pos[$i]$res_mut[$i] (+$res_mod[$i])";
		}
	}	
	else{
		$v = $res_mod[$i] - $res_prompt[$i];
		if($v < 0)	{
			
  			$e .= "<sup>$res_pos[$i]</sup>$res[$i]$v";
		}
		else{
  			$e .= "<sup>$res_pos[$i]</sup>$res[$i]+$v";
		}
		if($res_prompt[$i] != 0)	{
			$e .= ":$res_prompt[$i]";
		}
	}
	if($res[$i] =~ /m/i and abs($v - 16.0) < 0.5)	{
		$m_bMeOx = 1;
	}	
	if($res[$i] =~ /[st]/i and abs($v - 80.0) < 0.5)	{
		$m_bPhospho = 1;
	}	
	if(($res[$i] =~ /k/i or $i == 0) and abs($v - 43.0) < 0.5)	{
		$m_bUrea = 1;
	}	
}
my $colspan = 1;
if($bI)	{
	$colspan = 2;
}
if(defined $e)	{
  print "</tr><td><td colspan=\"3\">&nbsp;</td><td colspan=\"$colspan\" align=\"right\" valign=\"center\"><b>mods:</b></td><td align=\"left\" valign=\"top\" colspan=\"1\">";
  print "$e</td>";
}

print "</tr>";

if(defined $note_line)	{
	if($bI == 0)	{
		print "<tr><td colspan=\"6\">$note_line</td></tr>";
	}
	else	{
		print "<tr><td colspan=\"7\">$note_line</td></tr>";
	}
}
print "</table>";
print "<script type='text/javascript'>ajax_begin_pep_count_batch($divctr_start, $divctr);</script>\n";

if($bOmega)	{
#	print "<script type='text/javascript'>ajax_begin_omega_score($divctr_start, $divctr);</script>\n";
	print "<script type='text/javascript'>ajax_begin_omega_charge_score_batch($omegactr_start, $omegactr);</script>\n";
}

if($use_png == 0)	{
	draw_spectrum(\@mass_values,\@int_values,$pic_path,$mass);
	my $u = $cgi->url();
	$u =~ s/(http\:\/\/.+?)\/.+/$1/;
	print qq(
	<script language="javascript">
		writeSVG('$u$pic_path_url',300,800);
	</script>
	);
}
else	{
	draw_spectrum(\@mass_values,\@int_values,$pic_path,$mass);
	print "\n<img src=\"$pic_path_url\" height=\"300\" width=\"800\">\n";
}
stats_table();
ion_table();

if(length($annotation) > 3 or scalar(@fixed_mods) or scalar(@variable_mods))	{
	print qq(<br /><b>Residue modification sets tested:</b><br />\n);
	print qq(<table width="650pt" cellpadding="2" cellspacing="2" class="altb">\n);
	my $isCam = 0;
	if(scalar(@fixed_mods))	{
		print qq(<tr><td align="right" valign="top" width="150"><b>Complete mods:</b></td>
		<td align="left" valign="top"><ol style="margin-top: 0px;margin-bottom:0px;list-style-type: lower-roman;">);
		my $a = 1;
		foreach $_(@fixed_mods)	{
			$_ = markup_annotation($_);
			if(/Carbamidomethyl/)	{
				$isCam = 1;
			}
			print qq(<li>$_</li>\n);
			$a++;
		}
		print qq(</ol></td>\n</tr>\n);
	}
	if(scalar(@variable_mods))	{
		print qq(<tr><td align="right" valign="top" width="150"><b>Potential mods:</b></td>
		<td align="left" valign="top"><ol style="margin-top: 0px;margin-bottom:0px;list-style-type: lower-roman;">);
		my $a = 1;
		foreach $_(@variable_mods)	{
			$_ = markup_annotation($_);
			print qq(<li>$_</li>\n);
			$a++;
		}
		print qq(</ol></td>\n</tr>\n);
	}
	if(length($annotation) > 3)	{
		$annotation = markup_annotation($annotation);
		print qq(<tr><td align="right" valign="top" width="150"><b>Protein-specific PTMs:</b></td>
		<td align="left" valign="top"><ol style="margin-top: 0px;margin-bottom:0px;list-style-type: lower-roman;">
		<li>$annotation</li></ol></td>\n</tr>\n);
	}
	if($isCam)	{
		$annotation = qq(Ammonia-loss\@Q, Ammonia-loss\@C, Dehydrated\@E);
	}
	else	{
		$annotation = qq(Ammonia-loss\@Q, Dehydrated\@E);
	}
	print qq(<tr><td align="right" valign="top" width="150"><b>N-terminal:</b></td>
		<td align="left" valign="top">
		<ol style="margin-top: 0px;margin-bottom:0px;list-style-type: lower-roman;">
		<li> $annotation (peptide)</li>\n
		<li>ragged, Acetyl (protein)</li></ol></td>\n</tr>\n);
	print qq(</table>\n);
}
my @ks = keys(%accessions);
#LLAMA CHANGES
#print qq(<br /><b>Peptide found in the following );
#if(scalar(@ks) > 1)	{
#	print scalar(@ks);
#	print qq( proteins:</b>);
#}
#else	{
#	print qq(protein:</b>);
#}
#print qq(<div style="display:block" id="1001">);
#print qq(<table width="800" border="0" cellspacing="2" cellpadding="2">\n);
#my $alt = 1;
#my $pros = 1;
#foreach $_(@ks)	{
#	if($alt)	{
#		print qq(<tr class="alt"><td align="right" width="200" valign="top">$pros. 
#		<a href="/thegpm-cgi/protein.pl?path=$url&uid=$accessions{$_}&amp;label=$_&ltype=$ltype&homolog=$homolog&amp;proex=$proex" class="small_link">$_</a>
#		</td>);
#		$alt = 0;
#	}
#	else	{
#		print qq(<tr><td align="right" width="200" valign="top">$pros.
#		<a href="/thegpm-cgi/protein.pl?path=$url&uid=$accessions{$_}&amp;label=$_&ltype=$ltype&homolog=$homolog&amp;proex=$proex" class="small_link">$_</a>
#		</td>);
#		$alt = 1;
#	}
#	my $l = GetInfoNew($_);
#	print qq(<td align="left" valign="top">$l</td>);
#	$pros++;
#}
#print qq(</table>\n</div>\n);

output_gpm($gpm_show);
print_toggle_script();
if($error){
	print qq(
		<br /><br />
		<table><tr><td>
		The peptide information for spectra <b>$id</b>
		appears to be incomplete.</td></tr>
		<tr><td>Please <a href="mailto:contact\@beavisinformatics.ca?subject=peptide#$id">contact</a> us for help on resolving this issue.</td></tr>
		<tr><td>If possible, attach the file or portion of the file that caused this error message to appear.</td></tr>
		</table></div></div></body></html>
	);
}
print qq(<br><br><span class=\"small_label\"><a href="http://www.adobe.com/svg/viewer/install/main.html" target="_adobe">Adobe SVG plugin</a> may be required.
<BR><i>$file_version</i></span>);
my $svg_text;
open(SVG,"<..$pic_path");
while(<SVG>)	{
	$svg_text .= $_;
}
close(SVG);
my $m = 100.0*(2.0+int($m_fMassMax/100.0));
print qq(
<div style="display:none">
<form name="pepseek" method="GET" action="http://www.nwsr.manchester.ac.uk/cgi-bin/pepseeker/pepseek.pl" target="_info">
<input name="PEPTIDES" value="$domain">
<input name="which" value="AND">
<input name="ionscoresymbol" value="%3E%3D">
<input name="evaluesymbol" value="%3C">
<input name="pvaluesymbol" value="%3E">
<input type="submit" value="Search" name="login" >
</form>
</div>
<div style="display:none">
	<form action="http://$gvalue/thegpm-cgi/dblist_vpep.pl" method="post">
		<input type="hidden" name="seq" value="$domain">
		<input type="hidden" name="svg" value="$v">
		<input type="hidden" name="z" value="$charge">
		<input type="hidden" name="mvalue" value="$mass_line">
		<input type="hidden" name="ivalue" value="$int_line">
		<input type="hidden" name="svg_path" value="$pic_path">
		<input type="hidden" name="mass" value="$mass">
		<input type="hidden" name="m" value="$m">
		<textarea name="input_svg">$svg_text</textarea>
	</form>
</div>
);
my $mod_out;
$a = 0;
foreach $line(@res)	{
	$mod_out .= sprintf("%.3f:%.3f@%s[%i]\n",@res_mod[$a]+@res_mut[$a]-@res_prompt[$a],@res_prompt[$a],@res[$a],@res_pos[$a]);
	$a++;
}
my $form1;
my @ionkeys = keys(%m_IonTypes);
my $ion_types;
foreach $_(@ionkeys)	{
	if($m_IonTypes{$_})	{
		$ion_types .= $_;
	}
}

print qq(
<form action="/thegpm-cgi/peptide_studio.pl" method="POST" target="_studio">
<input type="hidden" name="seq" value="$domain">
<input type="hidden" name="mods" value="$mod_out">
<input type="hidden" name="fe" value="$m_fError">
<input type="hidden" name="fe_type" value="$m_pErrorType">
<input type="hidden" name="z" value="$charge">
<input type="hidden" name="start" value="$start">
<input type="hidden" name="end" value="$end">
<input type="hidden" name="mass" value="$mass">
<input type="hidden" name="mass_values" value="$mass_line">
<input type="hidden" name="int_values" value="$int_line">
<input type="hidden" name="path" value="$url">
<input type="hidden" name="label" value="$label">
<input type="hidden" name="uid" value="$uid_input">
<input type="hidden" name="homolog" value="$homolog">
<input type="hidden" name="id" value="$id">
<input type="hidden" name="$proex" value="proex">
<input type="hidden" name="ltype" value="$ltype">
<input type="hidden" name="ion_types" value="$ion_types">

</form>
);
print "$form1\n";
my $unmatched_text;
foreach $_(@spec_unmatched)	{
	$unmatched_text .= $_ . ",";
}

my $matched_text;
foreach $_(@spec_matched)	{
	$matched_text .= $_ . ",";
}
my $total_text;
foreach $_(@spec_total)	{
	$total_text .= $_ . ",";
}
my $mz_true = ($mass + $delta - $m_fProton)/$charge + $m_fProton;

print qq(
<div style="display:none">
	<form action="http://$gvalue/thegpm-cgi/peptide.pl" method="post">
		<textarea name="matched">$matched_text</textarea>
		<textarea name="unmatched">$unmatched_text</textarea>
		<textarea name="total">$total_text</textarea>
		<input type="hidden" value="unmatched" name="type">
		<input type="hidden" value="$mz_true" name="mz">
		<input type="hidden" value="$charge" name="z">
		<input type="hidden" value="$domain" name="seq">
		<input type="hidden" value="$path" name="path">
		<input type="hidden" value="$id" name="id">
	</form>
</div>
</div></div></body></HTML>
);


sub delta_table
{
	my ($width,$height,$xoffset,$yoffset) = @_;
	my $ppm = 0;
	if($m_pErrorType =~ /ppm/)	{
		$ppm = 1;
	}
	my $return = "";
	$return .= "<g id=\"delta\">\n";
	my $a = 1;
	my $error = $m_fError;
	my $line;
	my $top = $yoffset;
	my $bottom = $yoffset + $height;
	my $unit = $height/($m_lSequence - 1);
	my $center = $xoffset + $width/2.0;
	my $scale = $width/(4.0*$error);
	my $x = -1.0*$error*$scale + $center;
	my $y = $top + $unit/2;
	my $y1 = 0;
	my $y2 = 0;
	$return .= "<line class=\"delta\" x1=\"$x\" y1=\"$top\" x2=\"$x\" y2=\"$bottom\" />";
	$x += $error*$scale;
	$return .= "<line class=\"delta\" x1=\"$x\" y1=\"$top\" x2=\"$x\" y2=\"$bottom\" />";
	$x += $error*$scale;
	$return .= "<line class=\"delta\" x1=\"$x\" y1=\"$top\" x2=\"$x\" y2=\"$bottom\" />\n";
	my $sig_t = 0;
	my $sig_c = 0;
	my $dval;
	while($a < $m_lSequence)	{
		$error = 2*$m_fError;
		$y1 = $y - 2;
		$y2 = $y + 2;
		if($ppm)	{
			$error /= 1000000.0;
			$error *= @m_pfB[$a];
		}
		if(abs(@m_pfBd[$a]) < $error && $m_IonTypes{"b"})	{
			$x = -1*@m_pfBd[$a]*$scale + $center;
			$sig_c++;
			if($ppm)	{
				$x = -1*(1000000.0*@m_pfBd[$a]/@m_pfB[$a])*$scale + $center;
				$sig_t += (1000000.0*@m_pfBd[$a]/@m_pfB[$a])*(1000000.0*@m_pfBd[$a]/@m_pfB[$a]);				
			}
			else	{
				$sig_t += @m_pfBd[$a]*@m_pfBd[$a];
			}
			$dval = sprintf("%.03f b[%i] +1",@m_pfB[$a],$a);
			$line = sprintf("<line id=\"bline_$a\" xlink:title=\"$dval\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
			$return .= $line;
			$return .= "<title>$dval</title></line>\n";
		} 
		if($charge == 3 && $m_IonTypes{"b"}){
			if(abs(@m_pfBdx2[$a]) < $error)	{
					$x = -1*@m_pfBdx2[$a]*$scale + $center;
					$sig_c++;
					if($ppm)	{
						$x = -1*(1000000.0*@m_pfBdx2[$a]/(@m_pfB[$a]/2))*$scale + $center;
						$sig_t += (1000000.0*@m_pfBdx2[$a]/(@m_pfB[$a]/2))*(1000000.0*@m_pfBdx2[$a]/(@m_pfB[$a]/2));
					}
					else	{
						$sig_t += @m_pfBdx2[$a]*@m_pfBdx2[$a];
					}
					$dval = sprintf("%.03f b[%i] +2",$m_fProton+(@m_pfB[$a]-$m_fProton)/2,$a);
					$return .= sprintf("<line xlink:title=\"$dval\" id=\"blinez2_$a\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
					$return .= "<title>$dval</title></line>\n";
			} 
		}	
		elsif($charge == 4 && $m_IonTypes{"b"}){
			if(abs(@m_pfBdx2[$a]) < $error)	{
				$x = -1*@m_pfBdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfBdx2[$a]/(@m_pfB[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfBdx2[$a]/(@m_pfB[$a]/2))*(1000000.0*@m_pfBdx2[$a]/(@m_pfB[$a]/2));
				}
				else	{
					$sig_t += @m_pfBdx2[$a]*@m_pfBdx2[$a];
				}
				$dval = sprintf("%.03f b[%i] +2",$m_fProton+(@m_pfB[$a]-$m_fProton)/2,$a);
				$return .= sprintf("<line id=\"blinez2_$a\" xlink:title=\"$dval\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
			if(abs(@m_pfBdx3[$a]) < $error)	{
					$x = -1*@m_pfBdx3[$a]*$scale + $center;
					$sig_c++;
					if($ppm)	{
						$x = -1*(1000000.0*@m_pfBdx3[$a]/(@m_pfB[$a]/3))*$scale + $center;
						$sig_t += (1000000.0*@m_pfBdx3[$a]/(@m_pfB[$a]/3))*(1000000.0*@m_pfBdx3[$a]/(@m_pfB[$a]/3));
					}
					else	{
						$sig_t += @m_pfBdx3[$a]*@m_pfBdx3[$a];
					}
					$dval = sprintf("%.03f b[%i] +3",$m_fProton+(@m_pfB[$a]-$m_fProton)/2,$a);
					$return .= sprintf("<line id=\"blinez3_$a\" xlink:title=\"$dval\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
					$return .= "<title>$dval</title></line>\n";
			} 
		}	
		
		$error = 2*$m_fError;
		if($ppm)	{
			$error /= 1000000.0;
			$error *= @m_pfY[$a];
		}
		if(abs(@m_pfYd[$a]) < $error && $m_IonTypes{"y"})	{
			$x = -1*@m_pfYd[$a]*$scale + $center;
			$sig_c++;
			if($ppm)	{
				$x = -1*(1000000.0*@m_pfYd[$a]/@m_pfY[$a])*$scale + $center;
				$sig_t += (1000000.0*@m_pfYd[$a]/@m_pfY[$a])*(1000000.0*@m_pfYd[$a]/@m_pfY[$a]);
			}
			else	{
				$sig_t += @m_pfYd[$a]*@m_pfYd[$a];
			}
			$dval = sprintf("%.03f y[%i] +1",@m_pfY[$a],$m_lSequence-$a);
			$return .= sprintf("<line id=\"rline_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
			$return .= "<title>$dval</title></line>\n";
		} 
		if(($charge == 3 or $charge == 2) && $m_IonTypes{"y"}){
			if(abs(@m_pfYdx2[$a]) < $error)	{
				$x = -1*@m_pfYdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfYdx2[$a]/(@m_pfY[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfYdx2[$a]/(@m_pfY[$a]/2))*(1000000.0*@m_pfYdx2[$a]/(@m_pfY[$a]/2));
				}
				else	{
					$sig_t += @m_pfYdx2[$a]*@m_pfYdx2[$a];
				}
				$dval = sprintf("%.03f y[%i] +2",$m_fProton+(@m_pfY[$a]-$m_fProton)/2,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rlinez2_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}	
		if($charge == 4 && $m_IonTypes{"y"}){
			if(abs(@m_pfYdx2[$a]) < $error)	{
				$x = -1*@m_pfYdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfYdx2[$a]/(@m_pfY[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfYdx2[$a]/(@m_pfY[$a]/2))*(1000000.0*@m_pfYdx2[$a]/(@m_pfY[$a]/2));
				}
				else	{
					$sig_t += @m_pfYdx2[$a]*@m_pfYdx2[$a];
				}
				$dval = sprintf("%.03f y[%i] +2",$m_fProton+(@m_pfY[$a]-$m_fProton)/2,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rlinez2_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
			if(abs(@m_pfYdx3[$a]) < $error)	{
				$x = -1*@m_pfYdx3[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfYdx3[$a]/(@m_pfY[$a]/3))*$scale + $center;
					$sig_t += (1000000.0*@m_pfYdx3[$a]/(@m_pfY[$a]/3))*(1000000.0*@m_pfYdx3[$a]/(@m_pfY[$a]/3));
				}
				else	{
					$sig_t += @m_pfYdx3[$a]*@m_pfYdx3[$a];
				}
				$dval = sprintf("%.03f y[%i] +3",$m_fProton+(@m_pfY[$a]-$m_fProton)/3,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rlinez3_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}
			
		$error = 2*$m_fError;
		if($ppm)	{
			$error /= 1000000.0;
			$error *= @m_pfC[$a];
		}
		if(abs(@m_pfCd[$a]) < $error && $m_IonTypes{"c"})	{
			$x = -1*@m_pfCd[$a]*$scale + $center;
			$sig_c++;
			if($ppm)	{
				$x = -1*(1000000.0*@m_pfCd[$a]/@m_pfC[$a])*$scale + $center;
				$sig_t += (1000000.0*@m_pfCd[$a]/@m_pfC[$a])*(1000000.0*@m_pfCd[$a]/@m_pfC[$a]);
			}
			else	{
				$sig_t += @m_pfCd[$a]*@m_pfCd[$a];
			}
			$dval = sprintf("%.03f c[%i] +1",@m_pfC[$a],$a);
			$return .= sprintf("<line id=\"bline_$a\" xlink:title=\"$dval\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
			$return .= "<title>$dval</title></line>\n";
		} 
		if(($charge == 3 or $charge == 2) && $m_IonTypes{"c"}){
			if(abs(@m_pfCdx2[$a]) < $error)	{
				$x = -1*@m_pfCdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfCdx2[$a]/(@m_pfC[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfCdx2[$a]/(@m_pfC[$a]/2))*(1000000.0*@m_pfCdx2[$a]/(@m_pfC[$a]/2));
				}
				else	{
					$sig_t += @m_pfCdx2[$a]*@m_pfCdx2[$a];
				}
				$dval = sprintf("%.03f c[%i] +2",$m_fProton+(@m_pfC[$a]-$m_fProton)/2,$a);
				$return .= sprintf("<line id=\"blinez2_$a\" xlink:title=\"$dval\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}	
		if($charge == 4 && $m_IonTypes{"c"}){
			if(abs(@m_pfCdx2[$a]) < $error)	{
				$x = -1*@m_pfCdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfCdx2[$a]/(@m_pfC[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfCdx2[$a]/(@m_pfC[$a]/2))*(1000000.0*@m_pfCdx2[$a]/(@m_pfC[$a]/2));
				}
				else	{
					$sig_t += @m_pfCdx2[$a]*@m_pfCdx2[$a];
				}
				$dval = sprintf("%.03f c[%i] +2",$m_fProton+(@m_pfC[$a]-$m_fProton)/2,$a);
				$return .= sprintf("<line id=\"blinez2_$a\" xlink:title=\"$dval\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
			if(abs(@m_pfCdx3[$a]) < $error)	{
				$x = -1*@m_pfCdx3[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfCdx3[$a]/(@m_pfC[$a]/3))*$scale + $center;
					$sig_t += (1000000.0*@m_pfCdx3[$a]/(@m_pfC[$a]/3))*(1000000.0*@m_pfCdx3[$a]/(@m_pfC[$a]/3));
				}
				else	{
					$sig_t += @m_pfCdx3[$a]*@m_pfCdx3[$a];
				}
				$dval = sprintf("%.03f c[%i] +3",$m_fProton+(@m_pfC[$a]-$m_fProton)/3,$a);
				$return .= sprintf("<line id=\"blinez3_$a\" xlink:title=\"$dval\" class=\"bline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}	

		$error = 2*$m_fError;
		if($ppm)	{
			$error /= 1000000.0;
			$error *= @m_pfZ[$a];
		}
		if(abs(@m_pfZd[$a]) < $error && $m_IonTypes{"z"})	{
			$x = -1*@m_pfZd[$a]*$scale + $center;
			$sig_c++;
			if($ppm)	{
				$x = -1*(1000000.0*@m_pfZd[$a]/@m_pfZ[$a])*$scale + $center;
				$sig_t += (1000000.0*@m_pfZd[$a]/@m_pfZ[$a])*(1000000.0*@m_pfZd[$a]/@m_pfZ[$a]);
			}
			else	{
				$sig_t += @m_pfZd[$a]*@m_pfZd[$a];
			}
			$dval = sprintf("%.03f z[%i] +1",@m_pfZ[$a],$m_lSequence-$a);
			$return .= sprintf("<line id=\"rline_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
			$return .= "<title>$dval</title></line>\n";
		} 
		if(($charge == 3 or $charge == 2) && $m_IonTypes{"z"}){
			if(abs(@m_pfZdx2[$a]) < $error)	{
				$x = -1*@m_pfZdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfZdx2[$a]/(@m_pfZ[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfZdx2[$a]/(@m_pfZ[$a]/2))*(1000000.0*@m_pfZdx2[$a]/(@m_pfZ[$a]/2));
				}
				else	{
					$sig_t += @m_pfZdx2[$a]*@m_pfZdx2[$a];
				}
				$dval = sprintf("%.03f z[%i] +2",$m_fProton+(@m_pfZ[$a]-$m_fProton)/2,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rblinez2_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}	
		if($charge == 4 && $m_IonTypes{"z"}){
			if(abs(@m_pfZdx2[$a]) < $error)	{
				$x = -1*@m_pfZdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfZdx2[$a]/(@m_pfZ[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfZdx2[$a]/(@m_pfZ[$a]/2))*(1000000.0*@m_pfZdx2[$a]/(@m_pfZ[$a]/2));
				}
				else	{
					$sig_t += @m_pfZdx2[$a]*@m_pfZdx2[$a];
				}
				$dval = sprintf("%.03f z[%i] +2",$m_fProton+(@m_pfZ[$a]-$m_fProton)/2,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rlinez2_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
			if(abs(@m_pfZdx3[$a]) < $error)	{
				$x = -1*@m_pfZdx3[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfZdx3[$a]/(@m_pfZ[$a]/3))*$scale + $center;
					$sig_t += (1000000.0*@m_pfZdx3[$a]/(@m_pfZ[$a]/3))*(1000000.0*@m_pfZdx3[$a]/(@m_pfZ[$a]/3));
				}
				else	{
					$sig_t += @m_pfZdx3[$a]*@m_pfZdx3[$a];
				}
				$dval = sprintf("%.03f z[%i] +3",$m_fProton+(@m_pfZ[$a]-$m_fProton)/3,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rlinez3_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}	
		if(abs(@m_pfZ1d[$a]) < $error && $m_IonTypes{"z"})	{
			$x = -1*@m_pfZ1d[$a]*$scale + $center;
			$sig_c++;
			if($ppm)	{
				$x = -1*(1000000.0*@m_pfZ1d[$a]/@m_pfZ1[$a])*$scale + $center;
				$sig_t += (1000000.0*@m_pfZ1d[$a]/@m_pfZ1[$a])*(1000000.0*@m_pfZ1d[$a]/@m_pfZ1[$a]);
			}
			else	{
				$sig_t += @m_pfZ1d[$a]*@m_pfZ1d[$a];
			}
			$dval = sprintf("%.03f z[%i]+H +1",@m_pfZ1[$a],$m_lSequence-$a);
			$return .= sprintf("<line id=\"rline_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
			$return .= "<title>$dval</title></line>\n";
		} 
		if(($charge == 3 or $charge == 2) && $m_IonTypes{"z"}){
			if(abs(@m_pfZ1dx2[$a]) < $error)	{
				$x = -1*@m_pfZdx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfZ1dx2[$a]/(@m_pfZ1[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfZ1dx2[$a]/(@m_pfZ1[$a]/2))*(1000000.0*@m_pfZ1dx2[$a]/(@m_pfZ1[$a]/2));
				}
				else	{
					$sig_t += @m_pfZ1dx2[$a]*@m_pfZ1dx2[$a];
				}
				$dval = sprintf("%.03f z[%i]+H +2",$m_fProton+(@m_pfZ1[$a]-$m_fProton)/2,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rblinez2_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}	
		if($charge == 4 && $m_IonTypes{"z"}){
			if(abs(@m_pfZ1dx2[$a]) < $error)	{
				$x = -1*@m_pfZ1dx2[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfZ1dx2[$a]/(@m_pfZ1[$a]/2))*$scale + $center;
					$sig_t += (1000000.0*@m_pfZ1dx2[$a]/(@m_pfZ1[$a]/2))*(1000000.0*@m_pfZ1dx2[$a]/(@m_pfZ1[$a]/2));
				}
				else	{
					$sig_t += @m_pfZ1dx2[$a]*@m_pfZ1dx2[$a];
				}
				$dval = sprintf("%.03f z[%i]+H +2",$m_fProton+(@m_pfZ1[$a]-$m_fProton)/2,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rlinez2_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
			if(abs(@m_pfZ1dx3[$a]) < $error)	{
				$x = -1*@m_pfZ1dx3[$a]*$scale + $center;
				$sig_c++;
				if($ppm)	{
					$x = -1*(1000000.0*@m_pfZ1dx3[$a]/(@m_pfZ1[$a]/3))*$scale + $center;
					$sig_t += (1000000.0*@m_pfZ1dx3[$a]/(@m_pfZ1[$a]/3))*(1000000.0*@m_pfZ1dx3[$a]/(@m_pfZ1[$a]/3));
				}
				else	{
					$sig_t += @m_pfZdx3[$a]*@m_pfZ1dx3[$a];
				}
				$dval = sprintf("%.03f z[%i]+H +3",$m_fProton+(@m_pfZ1[$a]-$m_fProton)/3,$m_lSequence-$a);
				$return .= sprintf("<line id=\"rlinez3_$a\" xlink:title=\"$dval\" class=\"rline\" x1=\"%.0f\" y1=\"%.0f\" x2=\"%.0f\" y2=\"%.0f\" >",$x,$y1,$x,$y2);
				$return .= "<title>$dval</title></line>\n";
			} 
		}	

		$y += $unit;
		$a++;
	}
	$return .= "</g>";
	$return .= "<g id=\"tics\" class=\"tics\">\n";
	$x = $xoffset + 1;
	$y = $bottom + 5;
	$return .= "<line x1=\"$x\" y1=\"$bottom\" x2=\"$x\" y2=\"$y\" />\n";
	$x += $width/4 - 1;
	$return .= "<line x1=\"$x\" y1=\"$bottom\" x2=\"$x\" y2=\"$y\" />\n";
	$x += $width/4;
	$return .= "<line x1=\"$x\" y1=\"$bottom\" x2=\"$x\" y2=\"$y\" />\n";
	$x += $width/4;
	$return .= "<line x1=\"$x\" y1=\"$bottom\" x2=\"$x\" y2=\"$y\" />\n";
	$x = $width + $xoffset;
	$return .= "<line x1=\"$x\" y1=\"$bottom\" x2=\"$x\" y2=\"$y\" />\n";
	$return .= "<line x1=\"$xoffset\" y1=\"$bottom\" x2=\"$x\" y2=\"$bottom\" />\n";
	$return .= "</g>\n";
	$return .= "<g id=\"tic-labels\" class=\"label\">\n";
	$error = -1.0*$m_fError;
	$y = $bottom + 18;
	$x = $xoffset + $width/4 - 10;
	$return .= "<text x=\"$x\" y=\"$y\" font-family=\"verdana,arial\" font-size=\"10\">$error</text>\n";
	$error = 1.0*$m_fError;
	$x = $xoffset + $width - 10 - $width/4;
	$return .= "<text x=\"$x\" y=\"$y\" font-family=\"verdana,arial\" font-size=\"10\">$error</text>\n";
	$x = $center - 2.5*length("error (ppm)");
	$y += 15;
	if($ppm)	{
		$return .= "<text x=\"$x\" y=\"$y\" font-family=\"verdana,arial\" font-size=\"10\">error (ppm)</text>\n";
	}
	else	{
		$return .= "<text x=\"$x\" y=\"$y\" font-family=\"verdana,arial\" font-size=\"10\">error (Da)</text>\n";
	}
	if(not $sig_c)	{
		$sig_c = 1;
	}
	$m_sigma = sprintf("%.4f Da",sqrt($sig_t/$sig_c));
	if(sqrt($sig_t/$sig_c) > 0.01)	{
		$m_sigma = sprintf("%.2f Da",sqrt($sig_t/$sig_c));
	}
	if($ppm)	{
		$m_sigma = sprintf("%.0f ppm",sqrt($sig_t/$sig_c));
		if(sqrt($sig_t/$sig_c) < 10)	{
			$m_sigma = sprintf("%.1f ppm",sqrt($sig_t/$sig_c));
		}
	}
	$return .= "<!-- sigma=$m_sigma -->\n";
	$return .= "</g>\n";
	return $return;	
}

sub ion_table
{
	if(($m_IonTypes{"b"} == 1 or $m_IonTypes{"y"} == 1) and ($m_IonTypes{"c"} == 1 or $m_IonTypes{"z"} == 1))	{
		if($g_cid > $g_etd)	{
			ion_table_by();
			other_table();
		}
		elsif($g_etd > $g_cid)	{
			ion_table_cz();
		}
		else	{
			print qq(<table><tr><td>);
			ion_table_by();
			print qq(</td><td>);
			ion_table_cz();
			print qq(</td></tr></table>);
			other_table();
		}
		return;
	}
	if($m_IonTypes{"b"} == 1 or $m_IonTypes{"y"} == 1)	{
		ion_table_by();
	}
	if($m_IonTypes{"c"} == 1 or $m_IonTypes{"z"} == 1)	{
		ion_table_cz();
	}
	if($m_IonTypes{"b"} == 1 or $m_IonTypes{"y"} == 1)	{
		other_table();
	}
}

sub other_table
{
	print qq(<table width="400pt" cellpadding="1" cellspacing="1">
	<tr><td align="left" colspan="4"><b>Other observed ions:</b></td></tr>
	<tr>
		<td align="center" valign="middle"><b>Ion type</b></td>
		<td align="center" valign="middle"><b>m/z</b></td>
		<td align="center" valign="middle"><b>Ion type</b></td>
		<td align="center" valign="middle"><b>m/z</b></td>
	</tr>
);
	my @vs;
	my $v;
	my $length = scalar(@itable);
	my $c = 0;
	my $m;
	while($c < $length)	{
		print qq(<tr>\n);
		@vs = split /\t/,@itable[$c];
		$m = sprintf("%.3f",@vs[0]);
		print qq(<td align="center" valign="middle">@vs[1]</td><td align="center" valign="middle">$m</td>);
		$c++;
		if($c < $length)	{
			@vs = split /\t/,@itable[$c];
			$m = sprintf("%.3f",@vs[0]);
			print qq(<td align="center" valign="middle">@vs[1]</td><td align="center" valign="middle">$m</td>);
		}
		else	{
			print qq(<td align="center" valign="middle">&nbsp;</td><td align="center" valign="middle">&nbsp;</td>);
		}
		print "</tr>\n";
		$c++;
	}
	print qq(</table>\n);
}

sub ion_table_by
{
	if($charge > 2)	{
		print "<table width=\"700pt\" cellpadding=\"1\" cellspacing=\"1\">\n";
	}
	else	{
		print "<table width=\"375pt\" cellpadding=\"1\" cellspacing=\"1\">\n";
	}
	print "<tr>";
	
	print "<td align=\"center\"><b>bond</b></td>";
	print "<td align=\"center\" class=\"y0\"><sup>+1</sup>y</td>";
	print "<td align=\"center\" class=\"y17\"><sup>+1</sup>y<sup>-17</sup></td>";
	print "<td align=\"center\" class=\"y18\"><sup>+1</sup>y<sup>-18</sup></td>";
	print "<td align=\"center\" class=\"b0\"><sup>+1</sup>b</td>";
	print "<td align=\"center\" class=\"b17\"><sup>+1</sup>b<sup>-17</sup></td>";
	print "<td align=\"center\" class=\"b18\"><sup>+1</sup>b<sup>-18</sup></td>";
	if($charge > 2)	{
		print "<td align=\"center\" class=\"y0\"><sup>+2</sup>y</td>";
		print "<td align=\"center\" class=\"y17\"><sup>+2</sup>y<sup>-17</sup></td>";
		print "<td align=\"center\" class=\"y18\"><sup>+2</sup>y<sup>-18</sup></td>";
		print "<td align=\"center\" class=\"b0\"><sup>+2</sup>b</td>";
		print "<td align=\"center\" class=\"b17\"><sup>+2</sup>b<sup>-17</sup></td>";
		print "<td align=\"center\" class=\"b18\"><sup>+2</sup>b<sup>-18</sup></td>";
	}
	print "</tr>";
	my $value = 0;
	my $value1 = 0;
	my $a = 1;
	my $line;
	my $save;
	my $pvalue;
	my $p;
	while($a < $m_lSequence)	{
		print "<tr>";
		$line = @m_pSequence[$a-1];
		print "<td align=\"center\"><sup>$line</sup>$a</td>";
		$value = get_y($a);
		$save = $value;
		$line = sprintf("%.3f",$value);
		$pvalue = $m_fProton + ($value - $m_fProton)/2;
		if(is_found($value) and $charge == 2 and is_found($pvalue))	{
			print "<td align=\"center\" class=\"y0\"><sup>+1,2</sup>$line</td>\n";
		}
		elsif(is_found($value))	{
			print "<td align=\"center\" class=\"y0\">$line</td>\n";
		}
		elsif($charge == 2 and is_found($pvalue))	{
			print "<td align=\"center\" class=\"y0\"><sup>+2</sup>$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		$value -= $m_fAmmonia;
		$line = sprintf("%.3f",$value);
		if(is_found($value))	{
			print "<td align=\"center\" class=\"y17\">$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		$save -= $m_fWater;
		$line = sprintf("%.3f",$save);
		if(is_found($save))	{
			print "<td align=\"center\" class=\"y18\">$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		$value = get_b($a);
		$save = $value;
		$line = sprintf("%.3f",$value);
		$pvalue = $m_fProton + ($value - $m_fProton)/2;
		if(is_found($value) and $charge == 2 and is_found($pvalue))	{
			print "<td align=\"center\" class=\"b0\"><sup>+1,2</sup>$line</td>\n";
		}
		elsif(is_found($value))	{
			print "<td align=\"center\" class=\"b0\">$line</td>\n";
		}
		elsif($charge == 2 and is_found($pvalue))	{
			print "<td align=\"center\" class=\"b0\"><sup>+2</sup>$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		$value -= $m_fAmmonia;
		$line = sprintf("%.3f",$value);
		if(is_found($value))	{
			print "<td align=\"center\" class=\"b17\">$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		$save -= $m_fWater;
		$line = sprintf("%.3f",$save);
		if(is_found($save))	{
			print "<td align=\"center\" class=\"b18\">$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		if($charge > 2)	{
			$value = get_y($a);
			$save = $value;
			$value1 = $m_fProton + ($value-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"y0\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$value -= $m_fAmmonia;
			$value1 = $m_fProton + ($value-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"y17\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$save -= $m_fWater;
			$value1 = $m_fProton + ($save-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"y18\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$value = get_b($a);
			$save = $value;
			$value1 = $m_fProton + ($value-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"b0\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$value -= $m_fAmmonia;
			$value1 = $m_fProton + ($value-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"b17\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$save -= $m_fWater;
			$value1 = $m_fProton + ($save-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"b18\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
		}
		print "</tr>";
		$a++;
	}
	print "</table>";
}

sub ion_table_cz
{
	if($charge > 2)	{
		print "<table width=\"350pt\" cellpadding=\"1\" cellspacing=\"1\">\n";
	}
	else	{
		print "<table width=\"175pt\" cellpadding=\"1\" cellspacing=\"1\">\n";
	}
	print "<tr>";
	
	print "<td align=\"center\"><b>bond</b></td>";
	print "<td align=\"center\" class=\"y0\"><sup>+1</sup>z+1</td>";
	print "<td align=\"center\" class=\"y0\"><sup>+1</sup>z+2</td>";
	print "<td align=\"center\" class=\"b0\"><sup>+1</sup>c</td>";
	if($charge > 2)	{
		print "<td align=\"center\" class=\"y0\"><sup>+2</sup>z+1</td>";
		print "<td align=\"center\" class=\"y0\"><sup>+2</sup>z+2</td>";
		print "<td align=\"center\" class=\"b0\"><sup>+2</sup>c</td>";
	}
	print "</tr>";
	my $value = 0;
	my $value1 = 0;
	my $a = 1;
	my $line;
	my $save;
	my $pvalue;
	my $p1value;
	my $p;
	my $hydro = 1.007825035;
	while($a < $m_lSequence)	{
		print "<tr>";
		$line = @m_pSequence[$a-1];
		print "<td align=\"center\"><sup>$line</sup>$a</td>";
		$value = get_z($a);
		$p1value = $value + $hydro;
		$save = $value;
		$line = sprintf("%.3f",$value);
		$pvalue = $m_fProton + ($value - $m_fProton)/2;
		if(is_found($value))	{
			print "<td align=\"center\" class=\"y0\">$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		$line = sprintf("%.3f",$p1value);
		if(is_found($p1value))	{
			print "<td align=\"center\" class=\"y0\">$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		$value = get_c($a);
		$save = $value;
		$line = sprintf("%.3f",$value);
		$pvalue = $m_fProton + ($value - $m_fProton)/2;
		if(is_found($value) and $charge == 2 and is_found($pvalue))	{
			print "<td align=\"center\" class=\"b0\"><sup>+1,2</sup>$line</td>\n";
		}
		elsif(is_found($value))	{
			print "<td align=\"center\" class=\"b0\">$line</td>\n";
		}
		elsif($charge == 2 and is_found($pvalue))	{
			print "<td align=\"center\" class=\"b0\"><sup>+2</sup>$line</td>\n";
		}
		else	{
			print "<td align=\"center\">$line</td>\n";
		}
		if($charge > 2)	{
			$value = get_z($a);
			$save = $value;
			$value1 = $m_fProton + ($value-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"y0\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$value1 = $m_fProton + ($value+$hydro-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"y0\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$value = get_c($a);
			$save = $value;
			$value1 = $m_fProton + ($value-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
				print "<td align=\"center\" class=\"b0\">$line</td>\n";
			}
			else	{
				print "<td align=\"center\">$line</td>\n";
			}
			$value -= $m_fAmmonia;
			$value1 = $m_fProton + ($value-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
#				print "<td align=\"center\" class=\"b17\">$line</td>\n";
			}
			else	{
#				print "<td align=\"center\">$line</td>\n";
			}
			$save -= $m_fWater;
			$value1 = $m_fProton + ($save-$m_fProton)/2.0;
			$line = sprintf("%.3f",$value1);
			if(is_found($value1))	{
#				print "<td align=\"center\" class=\"b18\">$line</td>\n";
			}
			else	{
#				print "<td align=\"center\">$line</td>\n";
			}
		}
		print "</tr>";
		$a++;
	}
	print "</table>";
}


sub is_found
{
	my($mass) = @_;
	my $return = 0;
	my $a;
	my $value = 0;
	my $error = $m_fError;
	if($m_pErrorType =~ /ppm/)	{
		$error /= 1000000.0;
		$error *= $mass;
	}
	foreach $value(@mass_values)	{
		if(abs($mass - $value) < $error)	{
			$return = 1;
			last;
		}
	}
	return $return;
}

sub get_tics
{
	my ($max) = @_;
	$max = 100*(1 + int($max/1000));
	return $max;
}



sub draw_spectrum
{
	my($_m,$_i,$_p,$parent) = @_;
	my $tp = get_root() . $_p;
	if($use_png == 1)	{
		$tp = "temp.svg";
	}
	if(not open(OUTPUT,">$tp"))	{
		$tp = "../$_p";
		open(OUTPUT,">$tp") or die "$tp not found";
	}
print OUTPUT <<End_of_svg;
<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 20010904//EN"
"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="800" height="300">
<defs>
<style type="text/css">
<![CDATA[
	/* delta blips */
	line.bline{
		stroke:blue;
		stroke-width:4;
	}
	line.rline
	{
		stroke:red;
		stroke-width:4;
	}
	/* masses */
	text.mz
	{
		visibility:hidden;
		font-family:verdana,arial;
		font-size:10;
	}
	/* delta graph vertical lines*/
	line.delta
	{
		stroke:black;
		stroke-width:1;
	}
	/* general text */
	.label
	{
		font-family:verdana,arial; 
		font-size:10;
	}
	.tics
	{
		stroke:black; 
		stroke-width:2;
	}
	
	/*histogram peaks */
	.bluepk
	{
		stroke-width:2;
		stroke:blue;
	}
	.greenpk
	{
		stroke-width:2;
		stroke:green;
	}
	.aquapk
	{
		stroke-width:2;
		stroke:#00CCFF;
	}
	.purplepk
	{
		stroke-width:2;
		stroke:#CC66FF;
	}
	.blackpk
	{
		stroke-width:2;
		stroke:black;
		fill:none;
	}
	.redpk
	{
		stroke-width:2;
		stroke:red;
	}
	.orgpk
	{
		stroke-width:2;
		stroke:orange;
	}
	.yellowpk
	{
		stroke-width:2;
		stroke:#DDDD00;
	}
	.cyanpk
	{
		stroke-width:2;
		stroke:#999999;
	}
	/*fragment-o-gram peaks*/
	.bluefr
	{
		stroke-width:3;
		stroke:blue;
	}
	.greenfr
	{
		stroke-width:3;
		stroke:green;
	}
	.internalpk
	{
		stroke-width:2;
		stroke:#8b4513;
	}
	.aquafr
	{
		stroke-width:3;
		stroke:#00CCFF;
	}
	.purplefr
	{
		stroke-width:3;
		stroke:#CC66FF;
	}
	.redfr
	{
		stroke-width:3;
		stroke:red;
	}
	.orgfr
	{
		stroke-width:3;
		stroke:orange;
	}
]]></style>
</defs>

<g id="boundary" class="blackpk">
<rect x="50" y="55" width="550" height="200"/>
</g>	
End_of_svg
	set_ions();
	print OUTPUT draw_fragment($_m,$_i);
	my $a = 0;
	my $width = 550;
	my $height = 200;
	my $mass_width = 100.0*(2.0+int($m_fMassMax/100.0));
	my $tic_width = get_tics($mass_width);
	my $x1;
	my $y1;
	my $x2;
	my $y2;
	my $x_scale = $width/$mass_width;
	my $y_scale = 180.0/100.0;
	print OUTPUT "<g id=\"x-tics\" class=\"tics\">";
	while($a < $mass_width+1)	{
		$x1 = int(50+($a*$x_scale) + 0.5);
		$x2 = $x1;
		$y1 = 255;
		$y2 = 259;
		print OUTPUT "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\"/>\n";
		$a += $tic_width;
	}
	print OUTPUT "</g>\n";
	print OUTPUT "<g id=\"y-tics\" class=\"tics\">";
	$a = 0;
	while($a < 101)	{
		$x1 = 50;
		$x2 = 46;
		$y1 = 255-int(0.5+($a*$y_scale));
		$y2 = $y1;
		print OUTPUT "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\"/>\n";
		$a += 10;
	}
	print OUTPUT "</g>\n";
	print OUTPUT "<g id=\"x-tic-labels\" class=\"label\">";
	$a = 0;
	while($a < $mass_width+1)	{
		$x1 = int(50+($a*$x_scale) - length($a)*4 + 0.5);
		$y1 = 273;
		print OUTPUT "<text x=\"$x1\" y=\"$y1\" font-family=\"verdana,arial\" font-size=\"10\">$a</text>\n";
		$a += $tic_width;
	}
	print OUTPUT "</g>\n";
	print OUTPUT "<g id=\"y-tic-labels\" class=\"label\">";
	$a = 0;
	while($a < 101)	{
		$x1 = 45-length($a)*8;
		$y1 = int(259-$a*$y_scale + 0.5);
		print OUTPUT "<text x=\"$x1\" y=\"$y1\" font-family=\"verdana,arial\" font-size=\"10\">$a</text>\n";
		$a += 20;
	}
	print OUTPUT "</g>\n";
	
	print OUTPUT "<g id=\"x-axis-label\" class=\"label\">\n";
	print OUTPUT "<text  x=\"350\" y=\"290\" font-family=\"verdana,arial\" font-size=\"10\">m/z</text>\n";
	print OUTPUT "</g>\n";

	print OUTPUT "<g id=\"y-axis-label\" class=\"label\">\n";
	print OUTPUT "<text  x=\"6\" y=\"168\" font-family=\"verdana,arial\" font-size=\"10\">RI</text>\n";
	print OUTPUT "</g>\n";
	
	print OUTPUT "<g id=\"histogram\">\n";
	$a = 0;
	my $int_max = 1;
	foreach $line(@$_m)	{
		if(@$_i[$a] > $int_max)	{
			$int_max = @$_i[$a];
		}
		$a++;
	}
	$a = 0;
	$int_max = 100.0/$int_max;
	my $is_b = 0;
	my $is_y = 0;
	my $is_c = 0;
	my $is_z = 0;
	my $is_b17 = 0;
	my $is_y17 = 0;
	my $is_b18 = 0;
	my $is_y18 = 0;
	my $is_y_ox = 0;
	my $is_b_ox = 0;
	my $is_special = 0;
	my $is_internal = 0;
	my $is_a = 0;
	my $b = 0;
	my $error;
	
	my $c0;
	my $c17;
	my $c18;
	my $z0;
	my $z17;
	my $z18;

	my $b0;
	my $b17;
	my $b18;
	my $y0;
	my $y17;
	my $y18;
	my $red = "";
	my $orange = "";
	my $blue = "";
	my $green = "";
	my $aqua = "";
	my $purple = "";
	my $cyan = "";
	my $yellow = "";
	my @marked_i;
	my @marked_m;
	my @u_m;
	my @u_i;
	my $modOx = 63.99829;
	my $hydro = 1.007825035;
	my $modPhospho = 98.0;
	print OUTPUT "<g id=\"showhide\">\n";
	print OUTPUT "<text class=\"label\" x=\"45\" y=\"290\" font-family=\"verdana,arial\" font-size=\"10\">view:</text>\n";
	print OUTPUT "<polygon id=\"matched\" points=\"80,290 84,278 88,290\" style=\"fill:red\"/>\n"; 
	print OUTPUT "<polygon id=\"unmatched\" points=\"98,290 102,278 106,290\" style=\"fill:black\" />\n"; 
	print OUTPUT "<text class=\"mz\" id=\"mtext\" x=\"115\" y=\"290\" font-family=\"verdana,arial\" font-size=\"10\">matched\n";
	print OUTPUT "<set attributeName=\"visibility\" to=\"visible\" begin=\"matched.mouseover\" end=\"matched.mouseout\"/></text>\n";
	print OUTPUT "<text class=\"mz\" id=\"umtext\" x=\"115\" y=\"290\" font-family=\"verdana,arial\" font-size=\"10\">unmatched\n";
	print OUTPUT "<set attributeName=\"visibility\" to=\"visible\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></text>\n";
	print OUTPUT "</g>\n";
	foreach $line(@$_m)	{
		$is_y = 0;
		$is_b = 0;
		$is_y_ox = 0;
		$is_b_ox = 0;
		$is_c = 0;
		$is_z = 0;
		$is_y17 = 0;
		$is_b17 = 0;
		$is_y18 = 0;
		$is_b18 = 0;
		$is_a = 0;
		$is_special = 0;
		$b = 0;
		$error = $m_fError;
		if($m_pErrorType =~ /ppm/)	{
			$error /= 1000000.0;
			$error *= $line;
		}
		my $c;
		$is_internal = 0;
		my $d = 0;
		foreach $c(@m_pfI)	{
			if(abs($c - $line) < $error)	{
				$is_internal = 1;
				push(@itable,"$c\t" . @m_strI[$d]);
			}
			$d++;
		}
		my $lab;
		my $r;
		my $f;
		while($b < $m_lSequence && $is_a == 0 && $is_y == 0 && $is_y17 == 0 && $is_b == 0 && $is_b17 == 0 && $is_b18 == 0 && $is_y18 == 0)	{
			$y0 = @m_pfY[$b];
			$b0 = @m_pfB[$b];
			$c0 = @m_pfC[$b];
			$z0 = @m_pfZ[$b];
			$y17 = $y0 - $m_fAmmonia;
			$b17 = $b0 - $m_fAmmonia;
			$y18 = $y0 - $m_fWater;
			$b18 = $b0 -  $m_fWater;
			$f = $b;
			$r = $m_lSequence - $f + 1;
			if(abs($b0 - 28 - $line) < $error)	{
				$is_a = 1;
				my $p = $b;
				push(@itable,$b0 - 28 . "\ta[$p]");
				$lab = "a[$f] +1";
			}
			if(abs($y0 - $line) < $error)	{
				$is_y = 1;
				$lab = "y[$r] +1";
			}
			if($m_bMeOx and abs($y0 - $line - $modOx) < $error)	{
				$is_y_ox = 1;
				my $p = $m_lSequence - $b;
				push(@itable,$y0 - $modOx . "\ty[$p]-CH<sub>4</sub>SO");
				$lab = "y[$r]-CH4SO +1";
			}
			if($m_bMeOx and abs($b0 - $line - $modOx) < $error)	{
				$is_b_ox = 1;
				my $p = $b;
				$lab = "b[$f]-CH4SO +1";
				push(@itable,$b0 - $modOx . "\tb[$p]-CH<sub>4</sub>SO");
			}
			if($m_bPhospho and abs($y0 - $line - $modPhospho) < 1)	{
				$is_y_ox = 1;
				my $p = $m_lSequence - $b;
				$lab = "y[$r]-H3PO4 +1";
				push(@itable,$y0 - $modPhospho . "\ty[$p]-H<sub>3</sub>PO<sub>4</sub>");
			}
			if($m_bPhospho and abs($b0 - $line - $modPhospho) < 1)	{
				$is_b_ox = 1;
				my $p = $b;
				$lab = "b[$f]-H3PO4 +1";
				push(@itable,$b0 - $modPhospho . "\tb[$p]-H<sub>3</sub>PO<sub>4</sub>");
			}
			if(abs($c0 - $line) < $error)	{
				$is_c = 1;
				$lab = "c[$f] +1";
			}
			if(abs($z0 - $line) < $error)	{
				$lab = "z[$r] +1";
				$is_z = 1;
			}
			if(abs($z0+$hydro - $line) < $error)	{
				$is_z = 1;
				$lab = "z+1[$r] +1";
			}
			if(abs($y17 - $line) < $error)	{
				$is_y17 = 1;
				$lab = "y[$r]-NH3 +1";
			}
			if(abs($y18 - $line) < $error)	{
				$is_y18 = 1;
				$lab = "y[$r]-H2O +1";
			}
			if(abs($b0 - $line) < $error)	{
				$is_b = 1;
				$lab = "b[$f] +1";
			}
			if(abs($b17 - $line) < $error)	{
				$is_b17 = 1;
				$lab = "b[$f]-NH3 +1";
			}
			if(abs($b18 - $line) < $error)	{
				$is_b18 = 1;
				$lab = "b[$f]-H2O +1";
			}
			if($charge >= 2)	{
				if(abs($m_fProton + ($z0-$m_fProton)/2.0 - $line) < $error)	{
					$is_z = 1;
					$lab = "z[$r] +2";
				}
				if(abs($m_fProton + ($z0+$hydro-$m_fProton)/2.0 - $line) < $error)	{
					$is_z = 1;
					$lab = "z[$r]+H +2";
				}
				if(abs($m_fProton + ($c0-$m_fProton)/2.0 - $line) < $error)	{
					$is_c = 1;
					$lab = "c[$f] +2";
				}
				if(abs($m_fProton + ($y0-$m_fProton)/2.0 - $line) < $error)	{
					$is_y = 1;
					$lab = "y[$r] +2";
				}
				if(abs($m_fProton + ($y17-$m_fProton)/2.0 - $line) < $error)	{
					$is_y17 = 1;
					$lab = "y[$r]-NH3 +2";
				}
				if(abs($m_fProton + ($y18-$m_fProton)/2.0 - $line) < $error)	{
					$is_y18 = 1;
					$lab = "y[$r]-H2O +2";
				}
				if(abs($m_fProton + ($b0-$m_fProton)/2.0 - $line) < $error)	{
					$is_b = 1;
					$lab = "b[$f] +2";
				}
				if(abs($m_fProton + ($b17-$m_fProton)/2.0 - $line) < $error)	{
					$is_b17 = 1;
					$lab = "b[$f]-NH3 +2";
				}
				if(abs($m_fProton + ($b18-$m_fProton)/2.0 - $line) < $error)	{
					$is_b18 = 1;
					$lab = "b[$f]-H2O +2";
				}
				if($m_bMeOx and abs($m_fProton + ($y0 - $modOx - $m_fProton)/2.0 - $line) < $error)	{
					$is_y_ox = 1;
					$lab = "y[$r]-CH4SO +2";
					my $p = $m_lSequence - $b;
					push(@itable,$m_fProton + ($y0 - $modOx - $m_fProton)/2.0 . "\ty[$p]-CH<sub>4</sub>SO<sup>+2</sup>");
				}
				if($m_bMeOx and abs($m_fProton + ($b0 - $modOx - $m_fProton)/2.0 - $line) < $error)	{
					$is_b_ox = 1;
					$lab = "b[$f]-CH4SO +2";
					my $p = $b;
					push(@itable,$m_fProton + ($b0 - $modOx - $m_fProton)/2.0 . "\tb[$p]-CH<sub>4</sub>SO<sup>+2</sup>");
				}
				if($m_bPhospho and abs($m_fProton + ($y0 - $modPhospho - $m_fProton)/2.0 - $line) < $error)	{
					$is_y_ox = 1;
					my $p = $m_lSequence - $b;
					$lab = "y[$r]-H3PO4 +2";
					push(@itable,$m_fProton + ($y0 - $modPhospho - $m_fProton)/2.0 . "\ty[$p]-H<sub>3</sub>PO<sub>4</sub><sup>+2</sup>");
				}
				if($m_bPhospho and abs($m_fProton + ($b0 - $modPhospho - $m_fProton)/2.0 - $line) < $error)	{
					$is_b_ox = 1;
					my $p = $b;
					$lab = "b[$f]-H3PO4 +2";
					push(@itable,$m_fProton + ($b0 - $modPhospho - $m_fProton)/2.0 . "\tb[$p]-H<sub>3</sub>PO<sub>4</sub><sup>+2</sup>");
				} 
			}
			elsif($charge == 4)	{
				if(abs($m_fProton + ($c0-$m_fProton)/2.0 - $line) < $error)	{
					$is_c = 1;
					$lab = "c[$f] +2";
				}
				if(abs($m_fProton + ($z0-$m_fProton)/2.0 - $line) < $error)	{
					$is_z = 1;
					$lab = "z[$r] +2";
				}
				if(abs($m_fProton + ($z0+$hydro-$m_fProton)/2.0 - $line) < $error)	{
					$is_z = 1;
					$lab = "z[$r]+H +2";
				}
				if(abs($m_fProton + ($y0-$m_fProton)/2.0 - $line) < $error)	{
					$is_y = 1;
					$lab = "y[$r] +2";
				}
				if(abs($m_fProton + ($y17-$m_fProton)/2.0 - $line) < $error)	{
					$is_y17 = 1;
					$lab = "y[$r]-NH3 +2";
				}
				if(abs($m_fProton + ($y18-$m_fProton)/2.0 - $line) < $error)	{
					$is_y18 = 1;
					$lab = "y[$r]-H2O +2";
				}
				if(abs($m_fProton + ($b0-$m_fProton)/2.0 - $line) < $error)	{
					$is_b = 1;
					$lab = "b[$f] +2";
				}
				if(abs($m_fProton + ($b17-$m_fProton)/2.0 - $line) < $error)	{
					$is_b17 = 1;
					$lab = "b[$f]-NH3 +2";
				}
				if(abs($m_fProton + ($b18-$m_fProton)/2.0 - $line) < $error)	{
					$is_b18 = 1;
					$lab = "b[$f]-H2O +2";
				}
				if(abs($m_fProton + ($y0-$m_fProton)/3.0 - $line) < $error)	{
					$is_y = 1;
					$lab = "y[$r] +3";
				}
				if(abs($m_fProton + ($y17-$m_fProton)/3.0 - $line) < $error)	{
					$is_y17 = 1;
					$lab = "y[$r]-NH3 +3";
				}
				if(abs($m_fProton + ($y18-$m_fProton)/3.0 - $line) < $error)	{
					$is_y18 = 1;
					$lab = "y[$r]-H2O +3";
				}
				if(abs($m_fProton + ($b0-$m_fProton)/3.0 - $line) < $error)	{
					$is_b = 1;
					$lab = "b[$f] +3";
				}
				if(abs($m_fProton + ($b17-$m_fProton)/3.0 - $line) < $error)	{
					$is_b17 = 1;
					$lab = "b[$f]-NH3 +3";
				}
				if(abs($m_fProton + ($b18-$m_fProton)/3.0 - $line) < $error)	{
					$is_b18 = 1;
					$lab = "b[$f]-H2O +3";
				}
			}
			$b++;
		}
		$x1 = 50+int($line*$x_scale + 0.5);
		$x2 = $x1;
		$y1 = 255;
		$y2 = $y1-int(0.5+@$_i[$a]*$int_max*$y_scale);
		
		if($is_z == 1 && $m_IonTypes{"z"})	{
			$red .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"redpk\"><title>$line $lab</title>\n";
			$red .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		elsif($is_c == 1 && $m_IonTypes{"c"})	{
			$blue .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"bluepk\"><title>$line $lab</title>\n";
			$blue .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		
		elsif($is_y == 1 && $m_IonTypes{"y"})	{
			$red .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"redpk\"><title>$line $lab</title>\n";
			$red .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		elsif($is_b == 1 && $m_IonTypes{"b"})	{
			$blue .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"bluepk\"><title>$line $lab</title>\n";
			$blue .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		elsif($is_a == 1 && ($m_IonTypes{"b"} || $m_IonTypes{"y"}))	{
			$yellow .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"yellowpk\"><title>$line $lab</title>\n";
			$yellow .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		elsif($is_b18 == 1 && $m_IonTypes{"b"})	{
			$aqua .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"aquapk\"><title>$line $lab</title>\n";
			$aqua .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		elsif($is_y18 == 1 && $m_IonTypes{"y"})	{
			$purple .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"purplepk\"><title>$line $lab</title>\n";
			$purple .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		elsif($is_b17 == 1 && $m_IonTypes{"b"})	{
			$green .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"greenpk\"><title>$line $lab</title>\n";
			$green .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		elsif($is_y17 == 1 && $m_IonTypes{"y"})	{
			$orange .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"orgpk\"><title>$line $lab</title>\n";
			$orange .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
		
		elsif($is_b_ox && $m_IonTypes{"b"})	{
			$cyan .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"cyanpk\"><title>$line $lab</title>\n";
			$cyan .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
			$is_special=1;
		}
		elsif($is_y_ox && $m_IonTypes{"y"})	{
			$cyan .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"cyanpk\"><title>$line $lab</title>\n";
			$cyan .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
			$is_special=1;
		}
		elsif(is_special($line,$parent,$error,$charge))	{
			$cyan .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"cyanpk\"><title>$line $lab</title>\n";
			$cyan .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
			$is_special=1;
		}
		elsif($is_internal && $m_IonTypes{"b"} && $m_IonTypes{"y"})	{
			$cyan .= "<line xlink:title=\"$line $lab\" x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"internalpk\"><title>$line $lab</title>\n";
			$cyan .= "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/></line>\n";
		}
#		else	{
#			print OUTPUT "<line xlink:title=\"$line $lab\"x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"blackpk\"><title>$line $lab</title>\n";
#			print OUTPUT "<set attributeName=\"visibility\" to=\"hidden\" begin=\"matched.mouseover\" end=\"matched.mouseout\"/></line>\n";
#		}

		if($m_IonTypes{"y"} && ($is_y || $is_y17 || $is_y18) )	{
			push @marked_m,$line;
			push @marked_i,@$_i[$a];
			push @spec_matched,"$line @$_i[$a]";
			$spec_stats{'matched_ions'} = $spec_stats{'matched_ions'} + 1;
			$spec_stats{'matched_int'} = $spec_stats{'matched_int'} + @$_i[$a];
		}
		elsif($m_IonTypes{"b"} && ($is_b || $is_b17 || $is_b18) )	{
			push @marked_m,$line;
			push @marked_i,@$_i[$a];
			push @spec_matched,"$line @$_i[$a]";
			$spec_stats{'matched_ions'} = $spec_stats{'matched_ions'} + 1;
			$spec_stats{'matched_int'} = $spec_stats{'matched_int'} + @$_i[$a];
		}
		elsif($m_IonTypes{"c"} && ($is_c) )	{
			push @marked_m,$line;
			push @marked_i,@$_i[$a];
			push @spec_matched,"$line @$_i[$a]";
			$spec_stats{'matched_ions'} = $spec_stats{'matched_ions'} + 1;
			$spec_stats{'matched_int'} = $spec_stats{'matched_int'} + @$_i[$a];
		}
		elsif($m_IonTypes{"z"} && ($is_z) )	{
			push @marked_m,$line;
			push @marked_i,@$_i[$a];
			push @spec_matched,"$line @$_i[$a]";
			$spec_stats{'matched_ions'} = $spec_stats{'matched_ions'} + 1;
			$spec_stats{'matched_int'} = $spec_stats{'matched_int'} + @$_i[$a];
		}
		elsif($is_special || ($is_a && ($m_IonTypes{"b"} || $m_IonTypes{"y"})))	{
			push @marked_m,$line;
			push @marked_i,@$_i[$a];
			push @spec_matched,"$line @$_i[$a]";
			$spec_stats{'matched_ions'} = $spec_stats{'matched_ions'} + 1;
			$spec_stats{'matched_int'} = $spec_stats{'matched_int'} + @$_i[$a];
		}
		elsif($is_internal && ($m_IonTypes{"b"} && $m_IonTypes{"y"}) )	{
			push @marked_m,$line;
			push @marked_i,@$_i[$a];
			push @spec_matched,"$line @$_i[$a]";
			$spec_stats{'matched_ions'} = $spec_stats{'matched_ions'} + 1;
			$spec_stats{'matched_int'} = $spec_stats{'matched_int'} + @$_i[$a];
		}
		else	{
			push @u_m,$line;
			push @u_i,@$_i[$a];
			print OUTPUT "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"blackpk\">\n";
			print OUTPUT "<set attributeName=\"visibility\" to=\"hidden\" begin=\"matched.mouseover\" end=\"matched.mouseout\"/></line>\n";
			push @spec_unmatched,"$line @$_i[$a]";
			$spec_stats{'unmatched_ions'} = $spec_stats{'unmatched_ions'} + 1;
			$spec_stats{'unmatched_int'} = $spec_stats{'unmatched_int'} + @$_i[$a];
		}
		push @spec_total,"$line @$_i[$a]";
		$spec_stats{'total_ions'} = $spec_stats{'total_ions'} + 1;
		$spec_stats{'total_int'} = $spec_stats{'total_int'} + @$_i[$a];
		$a++;
	}
	print OUTPUT "$orange$green$blue$red$cyan$purple$aqua$yellow";
	print OUTPUT "</g>\n";
	
	print OUTPUT "<g id=\"masses\" class=\"label\">";
	$a = 0;
	my $e;
	my $length = scalar(@u_m);
	my $flip = 1;
	my $last_m = 0;
	my $last_i = 0;
	while($flip != 0)	{
		$a = 0;
		$flip = 0;
		while($a < $length - 1)	{
			if(@u_i[$a] < @u_i[$a+1])	{
				$last_m = @u_m[$a];
				$last_i = @u_i[$a];
				@u_m[$a] = @u_m[$a+1];
				@u_i[$a] = @u_i[$a+1];
				@u_m[$a+1] = $last_m;
				@u_i[$a+1] = $last_i;
				$flip++;
			}
			$a++;
		}		
	}
	my $i = 0;
	$a = 0;
	foreach $line(@marked_m)	{
		$e = sprintf("%.1f",$line);
		$x1 = 50+$line*$x_scale;
		$i = @marked_i[$a]*$int_max;
		$y1 = 253-($i*$y_scale);
		if($y1 > 100)	{
			$x1 += 3;
			print OUTPUT sprintf("<text transform=\"translate(%i,%i) rotate(-90)\" font-family=\"verdana,arial\" font-size=\"10\">$e\n",$x1,$y1);
		}
		else	{
			print OUTPUT sprintf("<text transform=\"translate($x1,$y1) rotate(0)\" font-family=\"verdana,arial\" font-size=\"10\">$e\n",$x1,$y1);
		}
		print OUTPUT "<set attributeName=\"visibility\" to=\"hidden\" begin=\"unmatched.mouseover\" end=\"unmatched.mouseout\"/>\n";
		print OUTPUT "</text>\n";
		$a++;
	}
	my $i = 0;
	$a = 0;
	foreach $line(@u_m)	{
		$e = sprintf("%.1f",$line);
		$x1 = 50+$line*$x_scale;
		$i = @u_i[$a]*$int_max;
		$y1 = 253-($i*$y_scale);
		if($y1 > 100)	{
			$x1 += 3;
			print OUTPUT "<text transform=\"translate($x1,$y1) rotate(-90)\" font-family=\"verdana,arial\" font-size=\"10\">$e\n";
		}
		else	{
			print OUTPUT "<text transform=\"translate($x1,$y1) rotate(0)\" font-family=\"verdana,arial\" font-size=\"10\">$e\n";
		}
		print OUTPUT "<set attributeName=\"visibility\" to=\"hidden\" begin=\"matched.mouseover\" end=\"matched.mouseout\"/>\n";
		print OUTPUT "</text>\n";
		
		$a++;
		if($a > 5)	{
			last;
		}
	}
	print OUTPUT "</g>\n";
	print OUTPUT delta_table(100,200,610,55);
	print OUTPUT "</svg>\n";
	close(OUTPUT);
	if($use_png == 1)	{
		$e = "java -jar c:\\xml-batik\\batik-rasterizer.jar -d " . get_root() . "$_p temp.svg >test.txt";
		system($e);
	}
}	

sub draw_fragment
{
	my($_m,$_i) = @_;
	my $output;
	my $a = 0;
	my $x1;
	my $y1;
	my $x2;
	my $y2;
	my $y_scale = 25.0/100.0;
	my $y_start;
	$output .= "<g id=\"fragment_sequence\" style=\"font-family:verdana,arial; font-size:12;\">\n";
	$x1 = 50;
	$y1 = 30;
	my $line;
	foreach $line(@m_pSequence)	{
		$output .= "<text x=\"$x1\" y=\"$y1\" font-family=\"verdana,arial\" font-size=\"12\">$line</text>\n";
		$x1 += 17;
	}
	$output .= "</g>\n";

	$output .=  "<g id=\"fragment_intensity\">\n";

	$a = 0;
	my $int_max = 1;
	my $mass_max = 0;
	foreach $line(@$_m)	{
		if(@$_i[$a] > $int_max)	{
			$int_max = @$_i[$a];
		}
		if($line > $mass_max)	{
			$mass_max = $line;
		}
		$a++;
	}
	$a = 0;
	$int_max = 100.0/$int_max;
	$m_fMassMax = $mass_max;
	my $is_type;
	my $is_value;
	my $b = 0;
	my $error;
	my $c0;
	my $z0;
	my $z10;
	my $b0;
	my $b17;
	my $b18;
	my $y0;
	my $y17;
	my $y18;
	$g_cid = 0;
	$g_etd = 0;
	my $hydro = 1.007825035;
	foreach $line(@$_m)	{
		## these vars ($is_y etc...) are now set to $b var if found and used for positioning
		## as well as the way to indicate whether or not that ion was found, 
		## instead of just setting $is_y etc... to 1 if found
		$is_type = "";
		$is_value = 0;		
		$b = 1;
		$error = $m_fError;
		if($m_pErrorType =~ /ppm/)	{
			$error /= 1000000.0;
			$error *= $line;
		}
		while($b < $m_lSequence){
			$y0 = @m_pfY[$b];
			$b0 = @m_pfB[$b];
			$c0 = @m_pfC[$b];
			$z0 = @m_pfZ[$b];
			$z10 = @m_pfZ1[$b];
			$y17 = $y0 - $m_fAmmonia;
			$b17 = $b0 - $m_fAmmonia;
			$y18 = $y0 - $m_fWater;
			$b18 = $b0 - $m_fWater;
			if(abs($y0 - $line) < abs(@m_pfYd[$b]) and $m_IonTypes{'y'})	{
				@m_pfYd[$b] = $y0 - $line;
			}
			if(abs($b0 - $line) < abs(@m_pfBd[$b]) and $m_IonTypes{'b'})	{
				@m_pfBd[$b] = $b0 - $line;
			}
			if(abs($c0 - $line) < abs(@m_pfCd[$b]) and $m_IonTypes{'c'})	{
				@m_pfCd[$b] = $c0 - $line;
			}
			if(abs($z0 - $line) < abs(@m_pfZd[$b]) and $m_IonTypes{'z'})	{
				@m_pfZd[$b] = $z0 - $line;
			}
			if(abs($z0 - $line) < abs(@m_pfZd[$b]) and $m_IonTypes{'z'})	{
				@m_pfZd[$b] = $z0 - $line;
			}
			if(abs($z10 - $line) < abs(@m_pfZ1d[$b]) and $m_IonTypes{'z'})	{
				@m_pfZ1d[$b] = $z10 - $line;
			}
			if(abs($c0 - $line) < $error and $m_IonTypes{'c'})	{
				$is_value = $b;
				$is_type = "c";
			}
			elsif(abs($z0 - $line) < $error and $m_IonTypes{'z'})	{
				$is_value = $b;
				$is_type = "z";
			}
			elsif(abs($z10 - $line) < $error and $m_IonTypes{'z'})	{
				$is_value = $b;
				$is_type = "z";
			}
			elsif(abs($y0 - $line) < $error and $m_IonTypes{'y'})	{
				$is_value = $b;
				$is_type = "y";
			}
			elsif(abs($y17 - $line) < $error and $m_IonTypes{'y'})	{
				$is_value = $b;
				$is_type = "y17";
			}
			elsif(abs($y18 - $line) < $error and $m_IonTypes{'y'})	{
				$is_value = $b;
				$is_type = "y18";
			}
			elsif(abs($b0 - $line) < $error and $m_IonTypes{'b'})	{
				$is_value = $b;
				$is_type = "b";
			}
			elsif(abs($b17 - $line) < $error and $m_IonTypes{'b'})	{
				$is_value = $b;
				$is_type = "b17";
			}
			elsif(abs($b18 - $line) < $error and $m_IonTypes{'b'})	{
				$is_value = $b;
				$is_type = "b18";
			}
			elsif($charge == 3  || ($charge == 2 and @m_pSequence[$b] == 'P'))	{
				if(abs($m_fProton + ($c0-$m_fProton)/2.0 - $line) < $error  and $m_IonTypes{'c'})	{
					if(abs($m_fProton + ($c0-$m_fProton)/2.0 - $line) < abs(@m_pfCdx2[$b]))	{
						@m_pfCdx2[$b] = $m_fProton + ($c0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "cx2";
				}
				if(abs($m_fProton + ($z0-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'z'})	{
					if(abs($m_fProton + ($z0-$m_fProton)/2.0 - $line) < abs(@m_pfZdx2[$b]))	{
						@m_pfZdx2[$b] = $m_fProton + ($z0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "z1x2";
				}
				if(abs($m_fProton + ($z10-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'z'})	{
					if(abs($m_fProton + ($z10-$m_fProton)/2.0 - $line) < abs(@m_pfZ1dx2[$b]))	{
						@m_pfZ1dx2[$b] = $m_fProton + ($z10-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "zx2";
				}
				if(abs($m_fProton + ($y0-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'y'})	{
					if(abs($m_fProton + ($y0-$m_fProton)/2.0 - $line) < abs(@m_pfYdx2[$b]))	{
						@m_pfYdx2[$b] = $m_fProton + ($y0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "yx2";
				}
				if(abs($m_fProton + ($y17-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'y'})	{
					$is_value = $b;
					$is_type = "y17x2";
				}
				if(abs($m_fProton + ($y18-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'y'})	{
					$is_value = $b;
					$is_type = "y18x2";
				}
				if(abs($m_fProton + ($b0-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'b'})	{
					if(abs($m_fProton + ($b0-$m_fProton)/2.0 - $line) < abs(@m_pfBdx2[$b]))	{
						@m_pfBdx2[$b] = $m_fProton + ($b0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "bx2";
				}
				if(abs($m_fProton + ($b17-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'b'})	{
					$is_value = $b;
					$is_type = "b17x2";
				}
				if(abs($m_fProton + ($b18-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'b'})	{
					$is_value = $b;
					$is_type = "b18x2";
				}
			}
			elsif($charge == 4)	{
				if(abs($m_fProton + ($c0-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'c'})	{
					if(abs($m_fProton + ($c0-$m_fProton)/2.0 - $line) < abs(@m_pfCdx2[$b]))	{
						@m_pfCdx2[$b] = $m_fProton + ($y0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "cx2";
				}
				if(abs($m_fProton + ($z0-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'z'})	{
					if(abs($m_fProton + ($z0-$m_fProton)/2.0 - $line) < abs(@m_pfZdx2[$b]))	{
						@m_pfZdx2[$b] = $m_fProton + ($z0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "zx2";
				}
				if(abs($m_fProton + ($z10-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'z'})	{
					if(abs($m_fProton + ($z10-$m_fProton)/2.0 - $line) < abs(@m_pfZ1dx2[$b]))	{
						@m_pfZ1dx2[$b] = $m_fProton + ($z10-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "z1x2";
				}
				if(abs($m_fProton + ($y0-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'y'})	{
					if(abs($m_fProton + ($y0-$m_fProton)/2.0 - $line) < abs(@m_pfYdx2[$b]))	{
						@m_pfYdx2[$b] = $m_fProton + ($y0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "yx2";
				}
				if(abs($m_fProton + ($y17-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'y'})	{
					$is_value = $b;
					$is_type = "y17x2";
				}
				if(abs($m_fProton + ($y18-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'y'})	{
					$is_value = $b;
					$is_type = "y18x2";
				}
				if(abs($m_fProton + ($b0-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'b'})	{
					if(abs($m_fProton + ($b0-$m_fProton)/2.0 - $line) < abs(@m_pfBdx2[$b]))	{
						@m_pfBdx2[$b] = $m_fProton + ($b0-$m_fProton)/2.0 - $line;
					}
					$is_value = $b;
					$is_type = "bx2";
				}
				if(abs($m_fProton + ($b17-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'b'})	{
					$is_value = $b;
					$is_type = "b17x2";
				}
				if(abs($m_fProton + ($b18-$m_fProton)/2.0 - $line) < $error and $m_IonTypes{'b'})	{
					$is_value = $b;
					$is_type = "b18x2";
				}

				if(abs($m_fProton + ($y0-$m_fProton)/3.0 - $line) < $error and $m_IonTypes{'y'})	{
					if(abs($m_fProton + ($y0-$m_fProton)/3.0 - $line) < abs(@m_pfYdx3[$b]))	{
						@m_pfYdx3[$b] = $m_fProton + ($y0-$m_fProton)/3.0 - $line;
					}
					$is_value = $b;
					$is_type = "yx3";
				}
				if(abs($m_fProton + ($y17-$m_fProton)/3.0 - $line) < $error and $m_IonTypes{'y'})	{
					$is_value = $b;
					$is_type = "y17x3";
				}
				if(abs($m_fProton + ($y18-$m_fProton)/3.0 - $line) < $error and $m_IonTypes{'y'})	{
					$is_value = $b;
					$is_type = "y18x3";
				}
				if(abs($m_fProton + ($b0-$m_fProton)/3.0 - $line) < $error and $m_IonTypes{'b'})	{
					if(abs($m_fProton + ($b0-$m_fProton)/3.0 - $line) < abs(@m_pfBdx3[$b]))	{
						@m_pfBdx3[$b] = $m_fProton + ($b0-$m_fProton)/3.0 - $line;
					}
					$is_value = $b;
					$is_type = "bx3";
				}
				if(abs($m_fProton + ($b17-$m_fProton)/3.0 - $line) < $error and $m_IonTypes{'b'})	{
					$is_value = $b;
					$is_type = "b17x3";
				}
				if(abs($m_fProton + ($b18-$m_fProton)/3.0 - $line) < $error and $m_IonTypes{'b'})	{
					$is_value = $b;
					$is_type = "b18x3";
				}
			}
			$b++;
		}
		$y1 = 25;
		$y_start = (@$_i[$a]*$int_max*$y_scale);
		if($y_start < 1.0)	{
			$y_start = 1.0;
		}
		if($is_type =~ /^[by]/)	{
			$g_cid++;
		}
		elsif($is_type =~ /^[cz]/)	{
			$g_etd++;
		}
		if($is_type eq "y18x3" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"purplefr\"/>\n";
		}
		if($is_type eq "y18x2" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"purplefr\"/>\n";
		}
		if($is_type eq "y18" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"purplefr\"/>\n";
		}
		if($is_type eq "b18x3" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"aquafr\"/>\n";
		}
		if($is_type eq "b18x2" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"aquafr\"/>\n";
		}
		if($is_type eq "b18" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"aquafr\"/>\n";
		}
		if($is_type eq "y17x3" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"orgfr\"/>\n";
		}
		if($is_type eq "y17x2" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"orgfr\"/>\n";
		}
		if($is_type eq "y17" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"orgfr\"/>\n";
		}
		if($is_type eq "b17x3" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"greenfr\"/>\n";
		}
		if($is_type eq "b17x2" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"greenfr\"/>\n";
		}
		if($is_type eq "b17" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"greenfr\"/>\n";
		}
		if($is_type eq "yx3" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"redfr\"/>\n";
		}
		if($is_type eq "yx2" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"redfr\"/>\n";
		}
		if($is_type eq "y" && $m_IonTypes{"y"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"redfr\"/>\n";
		}
		if($is_type eq "bx3" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"bluefr\"/>\n";
		}
		if($is_type eq "bx2" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"bluefr\"/>\n";
		}
		if($is_type eq "b" && $m_IonTypes{"b"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"bluefr\"/>\n";
		}
		if($is_type eq "c" && $m_IonTypes{"c"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"bluefr\"/>\n";
		}
		if($is_type eq "z" && $m_IonTypes{"z"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"redfr\"/>\n";
		}
		if($is_type eq "cx2" && $m_IonTypes{"c"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y_start + $y1;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"bluefr\"/>\n";
		}
		if($is_type eq "zx2" && $m_IonTypes{"z"})	{
			$x1 = 46+($is_value)*17;
			$x2 = $x1;
			$y2 = $y1 - $y_start;
			$output .= "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" class=\"redfr\"/>\n";
		}
		$a++;
	}
	$output .= "</g>\n";
	return $output;
}	


sub get_y
{
	my ($_l) = @_;
	my $value = $m_fY;
	if($_l > $m_lSequence)	{
		return 0.0;
	}
	$value += $m_fCleaveC - $m_fCleaveCdefault;
	my $a = $m_lSequence - 1;
	while($a >= $_l)	{
		$value += $m_pfAaMass->{@m_pSequence[$a]};
		$value += @m_pfMods[$a];
		$a--;
	}
	return $value + $m_fProton;
} 

sub get_b
{
	my ($_l) = @_;
	if($m_lSequence < $_l+1)	{
		return 0.0;
	}
	my $value = $m_fB;
	$value += $m_fCleaveN - $m_fCleaveNdefault;
	my $a = 0;
	while($a < $_l)	{
		$value += $m_pfAaMass->{@m_pSequence[$a]};
		$value += @m_pfMods[$a];
		$a++;
	}	
	return $value + $m_fProton;
} 

sub add_internals
{
	my $a = 1;
	my $ma;
	@m_pfI = ();
	@m_strI = ();
	my $p1;
	my $p2;
	while($a < $m_lSequence - 1)	{
		my $c = $m_lSequence - 1;
		while($c > $a)	{
			$b = $a;
			$ma = $m_fB;
			while($b <= $c)	{
				$ma += ($m_pfAaMass->{@m_pSequence[$b]} + @m_pfMods[$b]);
				$b++;
			}
			$ma += ($m_fY + $m_fProton - $m_fWater);
			$p1 = $a+1;
			$p2 = $c+1;
			$b -= 1;
			push(@m_pfI,$ma);
			push(@m_strI,"by[$p1-$p2]");
			$ma = $ma - $m_fB + $m_fA;
			$b -= 1;
			push(@m_pfI,$ma);
			push(@m_strI,"ay[$p1-$p2]");
			$c -= 1;
		}
		$a++;
	}
	$ma = ($mass - $m_fProton)/$charge + $m_fProton;
	push(@m_pfI,$ma);
	push(@m_strI,"parent<sup>+$charge</sup>");
	$ma = ($mass - $m_fProton-$m_fWater)/$charge + $m_fProton;
	push(@m_pfI,$ma);
	push(@m_strI,"(parent-H<sub>2</sub>O)<sup>+$charge</sup>");
	$ma = ($mass - $m_fProton-2.0*$m_fWater)/$charge + $m_fProton;
	push(@m_pfI,$ma);
	push(@m_strI,"(parent-2H<sub>2</sub>O)<sup>+$charge</sup>");
	$ma = ($mass - $m_fProton-$m_fAmmonia)/$charge + $m_fProton;
	push(@m_pfI,$ma);
	push(@m_strI,"(parent-NH<sub>3</sub>)<sup>+$charge</sup>");
	$ma = ($mass - $m_fProton-2.0*$m_fAmmonia)/$charge + $m_fProton;
	push(@m_pfI,$ma);
	push(@m_strI,"(parent-2NH<sub>3</sub>)<sup>+$charge</sup>");
	$a = scalar(@m_pSequence)-1;
	$ma = ($mass - $m_pfAaMass->{@m_pSequence[$a]} -  @m_pfMods[$a] - $m_fWater);
	$ma = ($ma - $m_fProton)/$charge + $m_fProton;
	push(@m_pfI,$ma);
	$p1 = @m_pSequence[$a];
	push(@m_strI,"(parent-$p1)<sup>+$charge</sup>");

}

sub get_z
{
	my ($_l) = @_;
	my $value = $m_fZ;
	if($_l > $m_lSequence)	{
		return 0.0;
	}
	$value += $m_fCleaveC - $m_fCleaveCdefault;
	my $a = $m_lSequence - 1;
	while($a >= $_l)	{
		$value += $m_pfAaMass->{@m_pSequence[$a]};
		$value += @m_pfMods[$a];
		$a--;
	}
	return $value + $m_fProton;
} 

sub get_c
{
	my ($_l) = @_;
	if($m_lSequence < $_l+1)	{
		return 0.0;
	}
	my $value = $m_fC;
	$value += $m_fCleaveN - $m_fCleaveNdefault;
	my $a = 0;
	while($a < $_l)	{
		$value += $m_pfAaMass->{@m_pSequence[$a]};
		$value += @m_pfMods[$a];
		$a++;
	}	
	return $value + $m_fProton;
} 

sub set_ions
{
	my $a = 0;
	@m_pfY = ();
	@m_pfB = ();
	@m_pfC = ();
	@m_pfZ = ();
	@m_pfZ1 = ();
	while($a < $m_lSequence)	{
		push @m_pfY,get_y($a);
		push @m_pfB,get_b($a);
		push @m_pfYd,1000.0;
		push @m_pfBd,1000.0;
		## +2 and +3 charged deltas
		push @m_pfYdx2,1000.0;
		push @m_pfBdx2,1000.0;
		push @m_pfYdx3,1000.0;
		push @m_pfBdx3,1000.0;

		push @m_pfZ,get_z($a);
		push @m_pfZ1,get_z($a)+1.007825035;
		push @m_pfC,get_c($a);
		push @m_pfZd,1000.0;
		push @m_pfZ1d,1000.0;
		push @m_pfCd,1000.0;
		## +2 and +3 charged deltas
		push @m_pfZdx2,1000.0;
		push @m_pfZ1dx2,1000.0;
		push @m_pfCdx2,1000.0;
		push @m_pfZdx3,1000.0;
		push @m_pfZ1dx3,1000.0;
		push @m_pfCdx3,1000.0;
		$a++;
	}
	add_internals();
}

sub set_mods
{
	my $a = 0;
	@m_pfY = ();
	@m_pfB = ();
	@m_pfMods = ();
	my $b = 0;
	my $l = scalar(@res);
	my $mod = 0;
	while($a < $m_lSequence)	{
		$b = 0;
		$mod = 0;
		while($b < $l)	{
			if((@res_pos[$b] - $start) == $a)	{
				$mod += @res_mod[$b];
			}
			$b++;
		}
		push @m_pfMods,$mod;
		push @m_pfY,get_y($a);
		push @m_pfB,get_b($a);	
		$a++;
	}
}

##rc - 20050225 - also check for phosphorylation
sub is_special
{
	my ($m,$mh,$e,$z) = @_;
	my $found=0;
	$mh = get_y(0);
	my $value = $mh - ($m_fWater*2.0 + 24.0);
	if(abs($m - $value) < $e)	{
		return 1;
	}
	$value = $mh - $m_fWater;
	if(abs($m - $value) < $e)	{
		push(@itable,$value . "\t(parent-H<sub>2</sub>O)");
		return 1;
	}
	$value = $mh - $m_fAmmonia;
	if(abs($m - $value) < $e)	{
		push(@itable,$value . "\t(parent-NH<sub>3</sub>O)");
		return 1;
	}
	$value = $m_fProton + (($mh - $m_fProton)/$z);
	if(abs($m - $value) < $e)	{
		push(@itable,$value . "\t(parent)<sup>+$z</sup>");
		return 1;
	}
	$value = $m_fProton + (($mh - 63.99829 - $m_fProton)/$z);
	if($m_bMeOx and abs($m - $value) < $e)	{
		push(@itable,$value . "\t(parent-CH<sub>4</sub>SO)<sup>+$z</sup>");
		return 1;
	}
	$value = $m_fProton + (($mh - 43.0 - $m_fProton)/$z);
	if($m_bUrea and abs($m - $value) < $e)	{
		push(@itable,$value . "\t(parent-(NH<sub>2</sub>)<sub>2</sub>CO)<sup>$z</sup>");
		return 1;
	}
	
	##check for phosphorylation (mod 80 +/- 1 da)
	my $mod_error=1;
	my $error=$m_fError;
	if($m_pErrorType =~ /ppm/)	{
		$error /= 1000000.0;
		$error *= $mh;
	}
	my $a=0;
	foreach (@res){
		if(/[STYR]/){
			if(abs($res_mod[$a] - 80 < $mod_error))	{
				$found=1;
				last;
			}
		}
		$a++;
	}
	if(!($found)){
		return 0;
	}
	$value = $m_fProton+(($mh-$m_fProton)-($res_mod[$a]+$m_fWater))/$charge;
	if(abs($value-$m)<$error){
		push(@itable,$value . "\t(parent-H<sub>3</sub>PO<sub>4</sub>)<sup>$charge</sup>");
		return 1;
	}
	$value = $m_fProton+(($mh-$m_fProton)-($res_mod[$a]+2.0*$m_fWater))/$charge;
	if(abs($value-$m)<$error){
		push(@itable,$value . "\t(parent-H<sub>5</sub>PO<sub>5</sub>)<sup>$charge</sup>");
		return 1;
	}
	
	if($charge > 1)	{
		return 0;
	}
	my $l = scalar(@seq) - 1;
	$value = $mh - $m_pfAaMass->{@seq[$l]};
	if(abs($m - $value) < $e)	{
		push(@itable,"$value\t(parent-@seq[$l])");
		return 1;
	}
	return 0;
}

sub generate_mgf
{
	my $matched = $cgi->param('matched');
	my $unmatched = $cgi->param('unmatched');
	my $total = $cgi->param('total');
	my $type = $cgi->param('type');
	my $mz = $cgi->param('mz');
	my $z = $cgi->param('z');
	my $seq = $cgi->param('seq');
	print qq(Content-type: text/plain\nContent-disposition: attachment; filename=$seq.mgf\n\n);
	print qq(SEARCH=MIS\r\nREPTYPE=Peptide\r\n## created by: $file_version\r\n##sequence: $seq\r\n);
	my @m = split /,/,$matched;
	my @u = split /,/,$unmatched;
	my @t = split /,/,$total;
	print "BEGIN IONS\r\n";
	print "PEPMASS=$mz\r\n";
	print "CHARGE=$z+\r\n";
	print "TITLE=matched ions, spectrum $spec_id, data set $spec_gpm, sequence $seq\r\n";
	foreach $_(@m)	{
		print "$_\r\n";
	}
	print "END IONS\r\n\r\n";
	print "BEGIN IONS\r\n";
	print "PEPMASS=$mz\r\n";
	print "CHARGE=$z+\r\n";
	print "TITLE=unmatched ions, spectrum $spec_id, data set $spec_gpm, sequence $seq\r\n";
	foreach $_(@u)	{
		print "$_\r\n";
	}
	print "END IONS\r\n\r\n";
	print "BEGIN IONS\r\n";
	print "PEPMASS=$mz\r\n";
	print "CHARGE=$z+\r\n";
	print "TITLE=all ions ($seq), spectrum $spec_id, data set $spec_gpm, sequence $seq\r\n";
	foreach $_(@t)	{
		print "$_\r\n";
	}
	print "END IONS\r\n\r\n";
}

sub stats_table
{
	print "<table cellpadding='2' cellspacing='2'>\n<tr>\n";
	my $rions = 100*$spec_stats{'matched_ions'}/$spec_stats{'total_ions'};
	my $rint = 100*$spec_stats{'matched_int'}/$spec_stats{'total_int'};
	print sprintf("<td><b>matched/total:</b></td><td width='100' align='center'># ions: %i%%</td><td width='100' align='center'>intensity: %i%%</td><td width='100' align='center'>&sigma;: %s</td></tr>\n",$rions,$rint,$m_sigma);
	print "</table>\n";
}

sub output_gpm
{
	my ($_s) = @_;
	#LLAMA CHANGES
	return;
	print "<a name=\"contributor\"></a><br /><table class=\"alt\" cellpadding=\"3\" cellspacing=\"3\" width=\"800\">";
	if(not $_s)	{
		print "<tr><td align=\"right\" width=\"100\" valign=\"top\"><b>Contributor:</b></td>\n";
		print "<td align=\"left\" valign=\"top\">";
		print "anonymous";
		print "</td>\n</tr></table>\n";
		return;
	}

	if($gpm_value{"name"})	{
		print "<tr><td align=\"right\" width=\"100\" valign=\"top\"><b>Contributor:</b></td>\n";
		print "<td align=\"left\" valign=\"top\">";
		my $v = decorate($gpm_value{"name"});
		$v =~ s/\<.+?\>//;
		print "$v";
		print "</td>\n</tr></table>\n";
	}


	print "<div id='named_attributes' style='display: block;'>\n";
	print "<table  class=\"alt\" border=\"0\" cellpadding=\"3\" cellspacing=\"3\" width=\"800\">\n";
	my @keys = keys(%gpm_value);
	my $k;
	my $v;
	foreach $k(@keys)	{
		if($k eq "name" or $k eq "add" or $k eq "anonymous" or $gpm_value{$k} =~ "none")	{
			next;
		}
		$v = $gpm_value{$k};
		$v = decorate($v);
		print "<tr><td align=\"right\" width=\"100\" valign=\"top\"><b>$k</b></td>\n";
		print "<td align=\"left\" valign=\"top\">";
		print $v;
		print "</td>\n</tr>\n";
	}
	print "</table>\n</div>";
}

sub print_toggle_script
{
	print qq(
	<script type="text/javascript">
	<!--
	function toggleBox(szDivID)
	{
		var obj = document.getElementById(szDivID);
		if(obj.style.display == "none"){
			obj.style.display = "block";
		}
		else{
			obj.style.display = "none";
		}
	}
	function hideBox(szDivID)
	{
		var obj = document.getElementById(szDivID);
		obj.style.display = "none";
	}
	function showBox(szDivID)
	{
		var obj = document.getElementById(szDivID);
		obj.style.display = "block";
	}
	// -->
	</script>
	);
}
