#!/usr/local/bin/perl
##
## Version 2004.02.03
## Version 2004.03.01
## Version 2004.06.08
## Version 2004.06.11
## Version 2004.07.19 - new style
## Version 2004.09.03
## 2004.10.19 - added get_cache_root and changed scd to sgd in info bar and protein report heading
## Version 2004.12.15 - updated regex in PostNcbi
## Version 2005.01.11 - added capability of accessing HUPO Nomenclature and UniProt from GetInfoBar 
## Version 2005.02.15 - added hash variable %g_ens for ensembl taxonomy to web directory translation used
##						in PostEnsembl and GetHrefEnsembl 
## Version 2005.02.22 - new function: set_aa($path) taken out of peptide,pgel,pchip and ptable.
##						it gets residue, H2O and NH3 mass values from output.xml if 
##						<group label="residue mass parameters" type="parameters"> exists in output.xml file.
##						It takes the output.xml file as a parameter and returns a hash ref containing the mass values.
## Version 2005.02.28 - new function: GetHrefGrid($label) adds capability of accesssing GRID site for yeast results in GetInfoBar.
## Version 2005.06.14 - added handling for UniProt id numbers in the description string
## Version 2005.07.07 - added handling for T. annulata genome sequences
## Version 2005.07.27 - altered ENSEMBL description handling for new ENSEMBL page format
## Version 2005.08.10 - corrected ENSEMBL description handling for new ENSEMBL page format
## Version 2005.08.11 - corrected ENSEMBL description handling problem of not recognizing a small number of results forms
## Version 2005.08.16 - added PRIDE link into the protein InfoBar
## Version 2005.09.07 - added the StartTop method
## Version 2005.11.03 - added support for protein composite validation display
## Version 2006.05.11 - added Fugu and Tetraodon ENSEMBL species
## Version 2006.05.17 - fix regex for IPI in GetInfoBar, add GetHrefPride to IPI and yeast in GetInfoBar
##					  - handle 'Archived Identifier' in PostEnsembl
##
## Version 2006.06.14 - javascript fix for svg requiring click to activate
## Version 2006.07.17 - added Monodelphis domestica & improved handling of T bruceii genome results
## Version 2006.07.18 - fixed problem with yeast Kegg links
## Version 2006.07.19 - updated link to GRID
## Version 2006.08.17 - added link to HMDB
## Version 2006.10.13 - added limited compatibility with PlasmoDB
## Version 2006.11.09 - added external references to obtain HGNC names
## Version 2007.02.16 - corrected C. elegans annotation change, necessary because of a change in the WormBase annotation system
## Version 2007.02.19 - added protein family information from ENSEMBL, if no description avaiable
## Version 2007.06.13 - corrected NCBI CGI URL
## Version 2007.07.22 - added 30 day check to update cached information requested by PostXxx
## Version 2007.09.14 - added code to follow ensembl archive links
## Version 2007.10.26 - updated rice TIGR posts and links
## Version 2007.11.27 - updated HGNC link
## Version 2008.04.11 - completed adding unlikely residue rules
## Version 2008.06.05 - added OMIM link to ENSMUSP pages
## version 2008.09.05 - modified cover_svg to fix a display bug in firefox that made clicking the svg graphic 
##	display protein.pl in the frame of the graphic.
## version 2009.05.20 - modified GetInfo to return MRM links for wormbase data
## version 2010.01.12 - modified cover_svg to return the percentage of the protein covered by peptides.
## version 2010.07.13 - modified getHrefMrm to call to the ajax_server_mrm.pl script to check 
##	the headers of a response for the presence of a label in the MRM Database instead of 
##	pulling back the whole response.
#
## common.pl
## Copyright (C) 2003-2005 Ronald C Beavis, all rights reserved
## The Global Proteome Machine 
## This software is a component of the X! proteomics software
## development project
##
## Use of this software governed by the Artistic license,
## as reproduced at http://www.opensource.org/licenses/artistic-license.php
##

## common.pl is used as an include file by many of the perl scripts in thegpm.
## it contains common functions used by these perl scripts. It has been merged with 
## get_info.pl functions


##  parameters: function dependant
##	called by: none
##	required by: check_spectra.pl, homolog.pl, model.pl, perform.pl, peptide.pl, plist.pl, protein.pl, transform_xml.pl
##	calls/links to: none

my %g_site = (
		"gpmdb_url"	=> "gpmdb.thegpm.org",
		"psyt_url"	=> "psyt.thegpm.org",
		"snap_url"	=> "snap.thegpm.org",
		"wiki_url"	=> "wiki.thegpm.org",
		"mrm_url"	=> "mrm.thegpm.org",
		"trace_url"	=> "ftp.thegpm.org",
		"h_channel"	=> "channel not set"
	);

my %g_trans = ("empty"	=> 1);
HurricaneSettings();

## global hash for ensembl taxonomy to web directory translation
## add an entry here when a new ensembl species is added
## $g_ens{ENSP} would return Homo_sapiens
my %g_meta=(		"AGAP"	=> 	"Anopheles_gambiae",
			"AT"	=>	"Arabidopsis_thaliana",
			"Bra"	=>	"Brassica_rapa",
			"ATCG"	=>	"Arabidopsis_thaliana",
			"ATMG"	=>	"Arabidopsis_thaliana",
			"LOC_Os"	=>	"Oryza_sativa",
			"LOC_Osp"	=>	"Oryza_sativa",
			"LOC_Osm"	=>	"Oryza_sativa",
			"FBpp" => "Drosophila_melanogaster",
			"DappuP" => "Daphnia_pulex",
			"PPA"  =>	"Pristionchus_pacificus",
			"BRADI"	=>	"Brachypodium_distachyon",
			"AAZ"	=>	"Trypanosoma_brucei",
			"AAQ"	=>	"Trypanosoma_brucei",
			"CAJ"	=>	"Trypanosoma_brucei",
			"EAN"	=>	"Trypanosoma_brucei",
			"CADAFUAP" =>	"Aspergillus_fumigatus"
	);
my %g_ens=(
		"ENSP" => "Homo_sapiens",
		"ENSSSCP" => "Sus_scrofa",
		"ENSMMUP" => "Macaca_mulatta",
		"ENSMODP" => "Monodelphis_domestica",
		"ENSCAFP" => "Canis_familiaris",
		"ENSMUSP" => "Mus_musculus",
		"ENSXETP" => "Xenopus_tropicalis",
		"ENSGALP" => "Gallus_gallus",
		"ENSRNOP" => "Rattus_norvegicus",
		"ENSANGP" => "Anopheles_gambiae",
		"AGAP" => "Anopheles_gambiae",
		"DappuP" => "Daphnia_pulex",
		"ENSDARP" => "Danio_rerio",
		"ENSBTAP" => "Bos_taurus",
		"ENSAPMP" => "Apis_mellifera",
		"ENSOCUP" => "Oryctolagus_cuniculus",
		"ENSCPOG" => "Cavia_porcellus",
		"ENSCPOP" => "Cavia_porcellus",
		"NEWSINFRUP" => "Fugu_rubripes",
		"ENSPTRP" => "Pan_troglodytes",
		"ENSFCAP" => "Felis_catus",
		"GSTENP"  => "Tetraodon_nigroviridis",
		"ENSTNIP" => "Tetraodon_nigroviridis",
		"ENSTRUP" => "Takifugu_rubripes",
		"ENSECAP" => "Equus_caballus",
		"ENSMEUP" => "Macropus_eugenii",
		"ENSCINP" => "Ciona_intestinalis",
		"FBpp" => "Drosophila_melanogaster",
		"ENSACAP" => "Anolis_carolinensis",
		"ENSTGUP" => "Taeniopygia_guttata",
		"ENSLAFP" => "Loxodonta_africana",
		"ENSMGAP" => "Meleagris_gallopavo",
		"ENSGACP"=> "Gasterosteus_aculeatus"
	);

my %g_hgnc;
my %g_pubmed;
my %g_swiss;
my %g_ipi;
my %g_info;
my %g_hgnc_types;
my $g_cookie_jar = $g_cookie_jar;

return 1;
use strict;
use CGI::Carp ('fatalsToBrowser');
use URI::Escape;
use HTTP::Cookies;

sub HurricaneSettings
{
	if(-e "hurricane.cfg")	{
		open(IN,"<hurricane.cfg");
		my @v;
		while(<IN>)	{
			chomp($_);
			@v = split /\t/,$_;
			if(@v[0] and @v[1])	{
				$g_site{@v[0]} = @v[1];
			}
		}
	}
}

sub GetGsite
{
	my ($v) = @_;
	return $g_site{$v};
}

sub DateCheck
{
	my ($_p,$_d) = @_;
	my ($atime,$mtime) = (stat($_p))[8,9];
	my $time = time;
	$atime = $time - $mtime;	
	my $update_period = 3600*24*$_d;
	if($atime > $update_period)	{
		return 0;
	}
	return 1;
}

sub set_aa
{
#	local $/ = "<\/group><\/group>";
	my $path=shift;
	my $value;
	my $m_pfAaMass;
	my @residue=();
	my @mass=();
	my $H2O=0;
	my $NH3=0;
	my $a=0;
	my $found=0;
	{
		$m_pfAaMass->{A} = 71.037110;
		$m_pfAaMass->{B} = 114.042930;
		$m_pfAaMass->{C} = 103.009190;
		$m_pfAaMass->{D} = 115.026940;
		$m_pfAaMass->{E} = 129.042590;
		$m_pfAaMass->{F} = 147.068410;
		$m_pfAaMass->{G} = 57.021460;
		$m_pfAaMass->{H} = 137.058910;
		$m_pfAaMass->{I} = 113.084060;
		$m_pfAaMass->{J} = 0.0;
		$m_pfAaMass->{K} = 128.094960;
		$m_pfAaMass->{L} = 113.084060;
		$m_pfAaMass->{M} = 131.040490;
		$m_pfAaMass->{N} = 114.042930;
		$m_pfAaMass->{O} = 0.0;
		$m_pfAaMass->{P} = 97.052760;
		$m_pfAaMass->{Q} = 128.058580;
		$m_pfAaMass->{R} = 156.101110;
		$m_pfAaMass->{S} = 87.032030;
		$m_pfAaMass->{T} = 101.047680;
		$m_pfAaMass->{U} = 150.953640;
		$m_pfAaMass->{V} = 99.068410;
		$m_pfAaMass->{W} = 186.079310;
		$m_pfAaMass->{X} = 111.060000;
		$m_pfAaMass->{Y} = 163.063330;
		$m_pfAaMass->{Z} = 128.058580;
		$m_pfAaMass->{H2O} = 18.01056470;
  		$m_pfAaMass->{NH3} = 17.02654911;
	}
	if(not open(INPUT,"<$path")){
		if(open(INPUT,"<$path.gz"))	{
			close(INPUT);
			system("gzip -d $path.gz");
			open(INPUT,"<$path");
		}
		else	{
			return $m_pfAaMass;
		}
	}

	while(<INPUT>){
  		if (/group label=\"residue mass parameters\" type=\"parameters\"/){ 
			my $group .= $_;
			$_ = <INPUT>;
			while($_ and not /\<\/group/)	{
				$group .= $_;
				$_ = <INPUT>;
			}
			$_ = $group;
  			$found=1;
  			(@residue) = /\<aa type=\"(\w)\" mass=\"\d+\.\d+\" \/>/g;
  			(@mass) = /\<aa type=\"\w\" mass=\"(\d+\.\d+)\" \/>/g;
  			($NH3) = /\<molecule type=\"NH3\" mass=\"(\d+\.\d+)\" \/>/g;
  			($H2O) = /\<molecule type=\"H2O\" mass=\"(\d+\.\d+)\" \/>/g;
  			foreach $value(@residue){
  				$m_pfAaMass->{$residue[$a]} = $mass[$a];
  				$a++;
  			}
  			$m_pfAaMass->{H2O} = $H2O;
  			$m_pfAaMass->{NH3} = $NH3;
  			close(INPUT);
  			last;
 
  		}
  	}
    return $m_pfAaMass;
}


## get_feature returns the attribute requested
## $line is: <group id="680" mh="1945.92" z="2" expect="4.3e-004" label="ENSP00000307705" type="model">
## to get the attribute for the charge, do this:
## $feature = "z";
## $charge = get_feature($line,$feature);
## returns 2

## updated as suggested by RH - 040719
sub get_feature
{
	my ($s,$f) = @_;
	my ($return) = $s =~ /$f=\"(.+?)\"/;
	if($return eq '0'){
		return 0;
	}
	return ($return || '_');
}


## PrintHeader prints the html header with the passed string used as the <TITLE>

## no changes - 040719
sub PrintHeader
{
my ($title) = @_;
print qq(Content-type: text/html


<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<HTML lang="en">
	<HEAD>
		<TITLE>GPM - $title</TITLE>
	<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />
	<meta http-equiv="X-UA-Compatible" content="chrome=1" />
      	<link type="text/css" rel="stylesheet" href="/llama-magic-html/tandem-style.css" />
      	<link type="text/css" rel="stylesheet" href="/llama-magic-html/tandem-style-print.css" media="print" />
	<link rel="SHORTCUT ICON" href="/favicon.ico" />

<script type="text/javascript" src="/llama-magic-html/writeSVG.js"></script>
);

my $server = get_server_name();
if($server =~ /thegpm\.org/)	{
	print qq(<script src="http://www.google-analytics.com/urchin.js" type="text/javascript">
	</script>
	<script type="text/javascript">
		_uacct = "UA-2678554-1";
		urchinTracker();
	</script>
	);
}

print qq(
<SCRIPT LANGUAGE="JavaScript" type="text/javascript"> 
<!--  
function toggleBox(szDivID)
{
	var obj = document.getElementById(szDivID);
	if(!obj)	{
		return;
	}
	if(obj.style.display == "none"){
		obj.style.display = "block";
	}
	else{
		obj.style.display = "none";
	}
}

function showBox(szDivID)
{
	var obj = document.getElementById(szDivID);
	if(!obj)	{
		return;
	}
	obj.style.display = "block";
}

function hideBox(szDivID)
{
	var obj = document.getElementById(szDivID);
	if(!obj)	{
		return;
	}
	obj.style.display = "none";
}

function NavTo(strGi)
{
	var strDbse = strGi;
	window.location.href=strDbse;  
}

//opens a popup to call compress script
function compressWindow(theURL,winName,features)
{
		var windowVariable = null;
	  	windowVariable = open(theURL,winName,features);
	  	return false;
}
function mover(message){window.status=message;return true;}
function mout(){window.status='';}

function hurrilinks()
{
	var l = 0;
	var size = document.links.length;
	while(l < size)	{
		var url = document.links[l].href;
		var myreg = /$g_site{"h_channel"}/;
		if(document.links[l].href.match(myreg))	{
			document.links[l].style.textDecoration = 'underline';
		}

		l++;
	}
}
//--> 
</SCRIPT>
</HEAD>
<BODY onLoad="hurrilinks();">
<div class="grad_gpm">
<div class="wrapper">
<span style="display:none">
<img src="/llama-magic-html/ffa_collapsed.gif" />
<img src="/llama-magic-html/ffa_expanded.gif" />
<img src="/llama-magic-html/waiting.gif" />
</span>
);
}

sub StartTop
{
	my ($file,$label) = @_;
	$file =~ s/\s+$//;
#LLAMA CHANGES
#	print qq(
#    		<table><tr><td align="center" valign="top">
#    		<a href="/index.html"><img src="/llama-magic-html/gpm.png" border="0"></a>
#    		</td><td>&nbsp;&nbsp;</td><td valign="middle" width="750">);
	print qq(
    		<table><tr><td align="center" valign="top">
		<a href="/index.html"><img src="/llama-magic-html/gpm.png" border="0"></a>
    		</td><td>&nbsp;&nbsp;</td><td valign="middle" width="750">);
#	if($label and $file)	{
#    		print qq(<i>$file: $label</i><BR><BR>);
#	}
#	elsif($label and not $file)	{
#   		 print qq(<i>$label</i><BR><BR>);
#	}
#	elsif(not $label and not $file)	{
#   		 print qq(<BR>);
#	}
#	else	{
#    		print qq(<i>$file</i><BR><BR>);
#	}
	print qq(<BR>);
}

## GetInfo recieves an accession number (ensembl, ncbi[gi], tigrat, mips etc..) and 
## opens the corresponding cache file and returns the html contained in that
## file. If it does not find a "cache.b" version of the cached file, it opens the full
## version, removes some html tags and then returns the revised html. Then it
## writes a new "cache.b" file.
## called by homolog.pl, protein.pl and plist.pl

## updated as suggested by RH - 040719

sub GetInfo	{
	my $val = shift;
  	my $cache = get_cache_root() . "/cache/";
  	my $html;
  	my $cpath;
  	my $ac;
  	my $l;
  	my $domains;
  	
	if($val =~ /\:reversed/)	{
		return "";
	}
 	if($val =~ /gb\|/)	{
	    $val =~ s/gi\|//;
	    $val =~ s/\|.*//g;
	    my $html = GetCache("$cache$val.ncbi.b");
	    return $html if ($html);
	    $html = GetCache("$cache$val.ncbi");
	    my @lines = split /\n/,$html;
	    my $a = 0;
	    my $s = scalar(@lines);
	    while(@lines[$a] !~ /^DEFINITION/ and $a < $s)	{
		$a++;
	    }
	    my $def;
	    my $extra = "";
	    my %doms;
	    if($a < $s)	{
	    	($def) = @lines[$a] =~ /^DEFINITION\s+?(.+)/;
		$a++;
		while($a < $s and @lines[$a] =~ /^\s/)	{
			($extra) = @lines[$a] =~ /^\s+(.+)/;
			$def .= " " . $extra;
			$a++;
		}
		while($a < $s and @lines[$a] !~ /^\s+\/gene\=\".+\"/)	{
			($extra) = @lines[$a] =~ /gene\=\"(.+?)\"/;
			$a++;
		}
		if($a < $s and @lines[$a] =~ /^\s+\/gene\=\".+\"/)	{
			($extra) = @lines[$a] =~ /gene\=\"(.+?)\"/;
			$def = "$extra, $def";
			$a++;
		}
	    }		
	    $a = 0;
		my $v;
		my $gene;
		my $br;
		my $href;
		$extra = "";
		while($a < $s)	{
			$v = @lines[$a];
			chomp($v);
			my $r;
			if($v =~ /\/region\_name\=/)	{
				($br) = $v =~ /\=\"*(.+)[\"\;]/;
				$gene = $br;
				$href = $br;
				$br = qq(<a href="http://www.ncbi.nlm.nih.gov/cdd?term=$href" target="_cdd" title="CDD domain $href">$br</a>);
				$r= "<br>$br";
				if($br)	{
					$a++;
					$v = @lines[$a];
					$v =~ s/\=\"/\=/;
					chomp($v);
					while($v !~ /\"/)	{
						$a++;
						$v .= " @lines[$a]";
						chomp($v);
					}
					if($v =~ /\/note\=/)	{
						($br) = $v =~ /\=(.+?)[\"]/;
						$br =~ s/\;.+//;
						my $len = 85 - length("$gene, ");
						$br =~ s/^(.{$len}).+/$1.../;
						$r .= ", $br";
						$doms{$r} = $doms{$r} + 1;
					}
				}
			}
			$a++;
		}
		my @ks = keys(%doms);
		if(scalar(@ks))	{
			foreach $extra(@ks)	{
				if($doms{$extra} > 1)	{
					$extra =~ s/\,/\($doms{$extra}\&times\;\)\,/;
				}
				$def .= $extra;
			}
		}	
	    WriteCache("$cache$val.ncbi.b",$def);
	    return ($def || "");
	}	
 	elsif($val =~ /gi\|/)	{
		$val =~ s/gi\|//;
		$val =~ s/\|.*//g;
	    my $html = GetCache("$cache$val.ncbi.b");
	    return $html if ($html);
	    $html = GetCache("$cache$val.ncbi");
	    my @lines = split /\n/,$html;
	    my $a = 0;
	    my $s = scalar(@lines);
	    while(@lines[$a] !~ /^DEFINITION/ and $a < $s)	{
		$a++;
	    }
	    my $def;
	    my $extra = "";
	    my %doms;
	    if($a < $s)	{
	    	($def) = @lines[$a] =~ /^DEFINITION\s+?(.+)/;
		$a++;
		while($a < $s and @lines[$a] =~ /^\s/)	{
			($extra) = @lines[$a] =~ /^\s+(.+)/;
			$def .= " " . $extra;
			$a++;
		}
		while($a < $s and @lines[$a] !~ /^\s+\/gene\=\".+\"/)	{
			($extra) = @lines[$a] =~ /gene\=\"(.+?)\"/;
			$a++;
		}
		if($a < $s and @lines[$a] =~ /^\s+\/gene\=\".+\"/)	{
			($extra) = @lines[$a] =~ /gene\=\"(.+?)\"/;
			$def = "$extra, $def";
			$a++;
		}
	    }		
	    $a = 0;
		my $v;
		my $br;
		my $href;
		my $gene;
		$extra = "";
		while($a < $s)	{
			$v = @lines[$a];
			chomp($v);
			my $r;
			if($v =~ /\/region\_name\=/)	{
				($br) = $v =~ /\=\"*(.+)[\"\;]/;
				$href = $br;
				$gene = $br;
				$br = qq(<a href="http://www.ncbi.nlm.nih.gov/cdd?term=$href" target="_cdd" title="CDD domain $href">$br</a>);
				$r= "<br>$br";
				if($br)	{
					$a++;
					$v = @lines[$a];
					$v =~ s/\=\"/\=/;
					chomp($v);
					while($v !~ /\"/)	{
						$a++;
						$v .= " @lines[$a]";
						chomp($v);
					}
					if($v =~ /\/note\=/)	{
						($br) = $v =~ /\=(.+?)[\"]/;
						$br =~ s/\;.+//;
						my $len = 85 - length("$gene, ");
						$br =~ s/^(.{$len}).+/$1.../;
						$r .= ", $br";
						
						$doms{$r} = $doms{$r} + 1;
					}
				}
			}
			$a++;
		}
		my @ks = keys(%doms);
		if(scalar(@ks))	{
			foreach $extra(@ks)	{
				if($doms{$extra} > 1)	{
					$extra =~ s/\,/ \<i\>\($doms{$extra}\&times\;\)\<\/i\>\,/;
				}
				$def .= $extra;
			}
		}	
	    WriteCache("$cache$val.ncbi.b",$def);
	    return ($def || "");
  	}	
 	elsif($val =~ /genedb\|/i)	{
		$val =~ s/genedb\|//i;
		$val =~ s/\|.*//g;
	    	my $html = GetCache("$cache$val.genedb.b");
	    	return $html if ($html);
	    	$html = GetCache("$cache$val.genedb");
		my $def = "no description available";
		if($html =~ /.+?\<th\>Product/si)	{
			$html =~ s/.+?\<th\>Product//si;
			($def) = $html =~ /\<span\>(.+?)\</si;
		}
		if(not scalar($def))	{
			$def = "no description available";
		}
	    	WriteCache("$cache$val.genedb.b",$def);
	    	return ($def || "");
  	}	
	elsif($val =~ /sp\|/)	{
	    $val =~ s/sp\|//;
	    $val =~ s/\|.*//g;
	    my $html = GetCache("$cache$val.sp.b");
	    $html =~ s/\w+\: \w+\=/ /g;
	    $html =~ s/Short\=/ /g;
	    $html =~ s/EC\=/EC /g;
	    $html = substr($html,0,256);
	    if(length($html) == 256)	{
		$html .= " ...";
	    }
	    return $html if ($html);
	    $html = GetCache("$cache$val.sp");
	    my @ls = split /\n/,$html;
	    my $z;
	    my $def;
	    foreach $z(@ls)	{
		if($z =~ /^DE /)	{
			chomp($z);
			$z =~ s/DE +//;
			$def .= "$z ";
		}
	    }
	    WriteCache("$cache$val.sp.b",$def);
	    $def =~ s/\w+\: \w+\=/ /g;
	    $def =~ s/Short\=/ /g;
	    $def =~ s/EC\=/EC /g;
	    $def = substr($def,0,256);
	    if(length($def) == 256)	{
		$def .= " ...";
	    }
	    return ($def || "");
  	}	
	elsif($val =~ /tr\|/)	{
	    $val =~ s/tr\|//;
	    $val =~ s/\|.*//g;
	    my $html = GetCache("$cache$val.sp.b");
	    $html =~ s/\w+\: \w+\=/ /g;
	    $html =~ s/Short\=/ /g;
	    $html =~ s/EC\=/EC /g;
	    $html = substr($html,0,256);
	    if(length($html) == 256)	{
		$html .= " ...";
	    }
	    return $html if ($html);
	    $html = GetCache("$cache$val.sp");
	    my @ls = split /\n/,$html;
	    my $z;
	    my $def;
	    foreach $z(@ls)	{
		if($z =~ /^DE /)	{
			chomp($z);
			$z =~ s/DE +//;
			$def .= "$z ";
		}
	    }
	    WriteCache("$cache$val.sp.b",$def);
	    $def =~ s/\w+\: \w+\=/ /g;
	    $def =~ s/Short\=/ /g;
	    $def =~ s/EC\=/EC /g;
	    $def = substr($def,0,256);
	    if(length($def) == 256)	{
		$def .= " ...";
	    }
	    return ($def || "");
  	}	
	elsif($val =~ /HIT[0-9]+/)	{
		$val =~ s/(HIT[0-9]+)\..+/$1/;
		my $html = GetCache("$cache$val.hit.b");
		return $html if($html);
		$html = GetCache("$cache$val.hit");
		my ($def) = $html =~ /Definition\<\/th\>.+?\<td\>(.+?)\<\/td\>/s;
		WriteCache("$cache$val.hit.b",$def);
	    return ($def || "");
	}	
	elsif($val =~ /tgo\|.+?\|/)	{
		$val =~ s/tgo\|(.+?)\|/$1/;
		my $html = GetCache("$cache$val.toxo.b");
		return $html if($html);
		$html = GetCache("$cache$val.toxo");
		my ($def) = $html =~ /\<h3\>.+?<\/h3\>\<h3\>(.+?)\<\/h3\>/s;
		WriteCache("$cache$val.toxo.b",$def);
	    return ($def || "");
	}	
	elsif($val =~ /IPI[0-9]+/ or /^IPI\:IPI[0-9]+/)	{
		$val =~ s/\.[0-9]*//;
		$val =~ s/^IPI\://;
		my $html = GetCache("$cache$val.ipi.b");
		return $html if($html);
		$html = GetCache("$cache$val.ipi");
		my ($def) = $html =~ /^DE\s+?(.+)$/m;
		WriteCache("$cache$val.ipi.b",$def);
	    return ($def || "");
	}	
	elsif($val =~ /tb[0-9]+/i)	{
		my $html = GetCache("$cache$val.tigrtb.b");
		return $html if($html);
		$html = GetCache("$cache$val.tigrtb");
		my ($def) = $html =~ /Gene Product Name.*?\<td.*?>(.*?)\<\/td>/si;
		WriteCache("$cache$val.tigrtb.b",$def);
	    return ($def || "");
	}	
	elsif($val =~ /TA[0-9][0-9][0-9][0-9][0-9]/)	{
		my $html = GetCache("$cache$val.genedbta.b");
		return $html if($html);
		$html = GetCache("$cache$val.genedbta");
		my ($def) = $html =~ /\<td\>\<b\>Product\<\/b\>\<\/td\>\s*\<td\>(.+?)\<\/td>/si;
	        $def =~ s/\(.*\)//gs;
		WriteCache("$cache$val.genedbta.b",$def);
	    return ($def || "");
	}	
	elsif($val =~ /osa1/)	{
		my $html = GetCache("$cache$val.tigrosa.b");
	        return $html if ($html);
		$html = GetCache("$cache$val.tigrosa");
		my ($def) = $html =~ /Gene Product Name.*?\<t[hd].*?\>(.+?)\<\/t[hd]\>/si;
		WriteCache("$cache$val.tirgosa.b",$def);
	        return ($def || '');
	}	
	elsif($val =~ /Y[A-Z]+?[0-9]+?[A-Z]/ or $val =~ /^[RQ]\d{4}/)	{
		my $html = GetCache("$cache$val.sgd.b");
		if(not $html =~ /SGD Protein Report\:/si)	{
			return $html if($html);
		}
		$html = GetCache("$cache$val.sgd");
		my $standard;
		if($html =~ /Standard name/si)	{
			$standard = $html;
			$standard =~ s/.*?Standard name.+?\<span class\=\'i\'\>(.*?)\<\/span\>.+/$1/si;
			$standard =~ s/\<.+?\>//sgi;
			$standard =~ s/\s*$//sgi;
		}
		else	{
			$standard = "";
		}
		my $description;
		if($html =~ /Description/s)	{
			$description = $html;
			$description =~ s/.*?Description\s*\<\/th\>.+?\<td.*?\>(.*?)\<\/td\>.+/$1/s;
			$description =~ s/\<.+?\>//sgi;
		}
		else	{
			$description = "";
		}
		if(length($standard) or length($description))	{
			$html = "$standard, $description";
		}
		else	{
			$html = "";
		}
		if(not $html =~ /SGD Protein Report\:/si)	{
			WriteCache("$cache$val.sgd.b",$html);
		}
		else	{
			$html = '';
		}
	    return ($html || '');
	}	
	elsif(0 and ($val =~ /^AGAP/ or ($val =~ /^FBpp/) or ($val =~ /^DappuP/) or $val =~ /^LOC\_Os/ or $val =~ /^AT[1-9C]G\d+/i  or isTb($val) or $val =~ /^Bra\d+\.\d/i
			or ($val =~ /[A-Z][0-9|A-Z]+?\.[0-9]/) or $val =~ /^BRADI/ or $val =~ /^SP[A-C][A-Z].*\-\d/ or $val =~ /^CADAFUAP/))  {
		my $html = GetCache("$cache$val.ensembl.b");
		$html =~ s/Source\:Uniprot\/SWISSPROT\;Acc\:(\w+)/\<a href=\"http\:\/\/www\.uniprot\.org\/uniprot\/$1\" target=\"_SPROT\"\>Source: UniProt $1\<\/a\>/s;
		if($html =~ /In Ensembl/)	{
			$html =~ s/In Ensembl .+?top of the page./no description available/si;
		}
		return $html if ($html);
		$html = GetCache("$cache$val.ensembl");
		my ($string) = $html =~ /\<p\>(.+?)\<\/p\>/si;
		if($string =~ /Description\: /si and not $string =~ /In Ensembl/si)	{
			$string =~ s/Description: //si;
		}
		else	{
			$string = "no description available";
		}
		my $ipr;
		if($html =~ /\<b\>Domains\<\/b/si)	{
			my ($v) = $html=~ /\<b\>Domains\<\/b\>(.+?)\<\/table/si;
			$v =~ s/\<td/\n\<td/gsi;
			$v =~ s/\<tr/\n\<tr/gsi;
			my @l = split /\n/,$v;
			my $td = 0;
			my $desc;
			my %iprs;
			my $i;
			my $count = 0;
			foreach $v(@l)	{
				if($v =~ /\<tr/)	{
					$td = 0;
					$desc = "";
					$count = 0;
				}
				elsif($v =~ /\<td/)	{
					$count++;
					if($count == 4)	{
						($desc) = $v =~ /\>(.+?)\</;
						$desc =~ s/\_/ /g;
					}
					elsif($count == 6)	{
						($i) = $v =~ /(\<a.+?\/a\>)/;
						if($i =~ /IPR/ and not $iprs{$i})	{
							$ipr .= "<br>$i $desc";
							$iprs{$i} = 1;
						}
						$count = 0;
					}
				}
			}
		}
		if(length($ipr))	{
			$string .="$ipr";
		}
		$string =~ s/\<br\>\&nbsp\;\&nbsp\;\&nbsp\;\&nbsp\;\[S/ [S/si;
		$string =~ s/Source\:UniProtKB\/Swiss\-Prot\;Acc\:(\w+)/\<a href=\"http\:\/\/www\.uniprot\.org\/uniprot\/$1\" target=\"_SPROT\"\>Source: UniProt $1\<\/a\>/si;
		WriteCache("$cache$val.ensembl.b",$string);
		return ($string || '');
	}
	elsif(($val =~ /ENS/) or ($val =~ /CG[0-9]+?-P[A-Z]/) or ($val =~ /NEWSINFRUP/) 
		or ($val =~ /GSTENP/) or $val =~ /^AGAP/ or ($val =~ /^FBpp/) or ($val =~ /^DappuP/) or $val =~ /^LOC\_Os/ or $val =~ /^AT[1-9C]G\d+/i or $val =~ /^Bra\d+\.\d/i 
		or isTb($val) 
			or ($val =~ /[A-Z][0-9|A-Z]+?\.[0-9]/) or $val =~ /^BRADI/ or $val =~ /^SP[A-C][A-Z].*\-\d/ or $val =~ /^CADAFUAP/)  {
		my $html = GetCache("$cache$val.ensembl.b");
		$html =~ s/Source\:Uniprot\/SWISSPROT\;Acc\:(\w+)/\<a href=\"http\:\/\/www\.uniprot\.org\/uniprot\/$1\" target=\"_SPROT\"\>Source: UniProt $1\<\/a\>/s;
		$html =~ s/Source\:HGNC Symbol\;Acc\:(\d+)/\<a href=\"http\:\/\/www\.genenames\.org\/data\/hgnc_data\.php\?hgnc_id\=$1\" target=\"_SPROT\"\>Source: HGNC $1\<\/a\>/s;
		if($html =~ /In Ensembl/)	{
			$html =~ s/In Ensembl .+?top of the page./no description available/si;
		}
		return $html if ($html);
		$html = GetCache("$cache$val.ensembl");
		if($html =~ /\<\!\-ENSEMBL51\-\-\>/)	{
			my ($string) = $html =~ /\<p\>(.+?)\<\/p\>/si;
			if($string =~ /Description\: /si and not $string =~ /In Ensembl/si)	{
				$string =~ s/Description: //si;
			}
			else	{
				$string = "no description available";
			}
			my $ipr;
			if($html =~ /\<b\>Domains\<\/b/si)	{
				my ($v) = $html=~ /\<b\>Domains\<\/b\>(.+?)\<\/table/si;
				$v =~ s/\<td/\n\<td/gsi;
				$v =~ s/\<tr/\n\<tr/gsi;
				my @l = split /\n/,$v;
				my $td = 0;
				my $desc;
				my %iprs;
				my $i;
				my $count = 0;
				foreach $v(@l)	{
					if($v =~ /\<tr/)	{
						$td = 0;
						$desc = "";
						$count = 0;
					}
					elsif($v =~ /\<td/)	{
						$count++;
						if($count == 4)	{
							($desc) = $v =~ /\>(.+?)\</;
							$desc =~ s/\_/ /g;
						}
						elsif($count == 6)	{
							($i) = $v =~ /(\<a.+?\/a\>)/;
							if($i =~ /IPR/ and not $iprs{$i})	{
								$ipr .= "<br>$i $desc";
								$iprs{$i} = 1;
							}
							$count = 0;
						}
					}
				}
			}
			if(length($ipr))	{
				$string .="$ipr";
			}
			$string =~ s/\<br\>\&nbsp\;\&nbsp\;\&nbsp\;\&nbsp\;\[S/ [S/si;
			$string =~ s/Source\:UniProtKB\/Swiss\-Prot\;Acc\:(\w+)/\<a href=\"http\:\/\/www\.uniprot\.org\/uniprot\/$1\" target=\"_SPROT\"\>Source: UniProt $1\<\/a\>/si;
			$string =~ s/Source\:HGNC Symbol\;Acc\:(\d+)/\<a href=\"http\:\/\/www\.genenames\.org\/data\/hgnc_data\.php\?hgnc_id\=$1\" target=\"_SPROT\"\>Source: HGNC $1\<\/a\>/si;
			WriteCache("$cache$val.ensembl.b",$string);
			return ($string || '');
		}
		elsif($html =~ /\<\!\-ENSEMBL60\-\-\>/)	{
			my ($string) = $html =~ /\<dt\>(.+?)\<\/dt\>/si;
			if($string =~ /Description/si and not $string =~ /In Ensembl/si)	{
				($string) = $html =~ /\<dd\>(.*?)\<\/dd\>/si;
				if(not scalar($string))	{
					$string = "no protein text annotation available";
				}
			}
			else	{
				$string = "no description available";
			}
			my $ipr;
			if($html =~ /\<b\>Domains\<\/b/si)	{
				my ($v) = $html=~ /\<b\>Domains\<\/b\>(.+?)\<\/table/si;
				$v =~ s/\<td/\n\<td/gsi;
				$v =~ s/\<tr/\n\<tr/gsi;
				my @l = split /\n/,$v;
				my $td = 0;
				my $desc;
				my %iprs;
				my $i;
				my $count = 0;
				foreach $v(@l)	{
					if($v =~ /\<tr/)	{
						$td = 0;
						$desc = "";
						$count = 0;
					}
					elsif($v =~ /\<td/)	{
						$count++;
						if($count == 4)	{
							($desc) = $v =~ /\>(.+?)\</;
							$desc =~ s/\_/ /g;
						}
						elsif($count == 6)	{
							($i) = $v =~ /(\<a.+?\/a\>)/;
							if($i =~ /IPR/ and not $iprs{$i})	{
								$ipr .= "<br>$i $desc";
								$iprs{$i} = 1;
							}
							$count = 0;
						}
					}
				}
			}
			if(length($ipr))	{
				$string .="$ipr";
			}
			$string =~ s/\<br\>\&nbsp\;\&nbsp\;\&nbsp\;\&nbsp\;\[S/ [S/si;
			$string =~ s/Source\:UniProtKB\/Swiss\-Prot\;Acc\:(\w+)/\<a href=\"http\:\/\/www\.uniprot\.org\/uniprot\/$1\" target=\"_SPROT\"\>Source: UniProt $1\<\/a\>/si;
			$string =~ s/Source\:HGNC Symbol\;Acc\:(\d+)/\<a href=\"http\:\/\/www\.genenames\.org\/data\/hgnc_data\.php\?hgnc_id\=$1\" target=\"_SPROT\"\>Source: HGNC $1\<\/a\>/si;
			WriteCache("$cache$val.ensembl.b",$string);
			return ($string || '');
		}
		else	{
			my ($description) = $html =~ /Description<\/th>.+?>(.+?)<br/si;
			my $domains;
			if(not length($description))	{
				($description) = $html =~ /Description.+?\<p\>(.*?)\<\/p/si;
			}
			else	{
				$description =~ s/<small>(.+?)<\/small>/$1/;
			}
			if(not length($description))	{
				($description) = $html =~ /Protein Family.+?\<p\>.*:(.*?)\<br.*?\<\/p/si;
				if(length($description))	{
					$description = "No description available.<br />Protein Family: " . $description;
				}
			}
			$description =~ s/Source\:Uniprot\/SWISSPROT\;Acc\:(\w+)/\<a href=\"http\:\/\/www\.uniprot\.org\/uniprot\/$1\" target=\"_SPROT\"\>Source: UniProt $1\<\/a\>/s;
			$description =~ s/Source\:HGNC Symbol\;Acc\:(\d+)/\<a href=\"http\:\/\/www\.genenames\.org\/data\/hgnc_data\.php\?hgnc_id\=$1\" target=\"_SPROT\"\>Source: HGNC $1\<\/a\>/s;
			my ($domainstring) = $html =~ /InterPro<\/th>(.+?)<\/table>/si;
			if(not length($domainstring))	{
				($domainstring) = $html =~ /InterPro.+?<\/th>(.+?)<\/table>/si;
				$domainstring =~ s/\<\/td\>/&nbsp;/gsi;
				$domainstring =~ s/\<td\>/&nbsp;/gsi;
				my @domains = $domainstring =~ m|(<a href=\"http\://www\.ebi\.ac\.uk\/interpro\/.+?)- |gsi;
				$domains = join "<br>&nbsp\;\n", @domains;
				$domains =~ s/<([c-zC-Z]).*?>//gs;
				$domains =~ s/<\/([c-zC-Z]).*?>//gs;
			}
			else	{
				$domainstring =~ s/<\/td>/&nbsp;/gsi;
				my @domains = $domainstring =~ m|(<A\sHREF=\"http\://www\.ebi\.ac\.uk\/interpro\/.+?)- |gsi;
				$domains = join "<br>&nbsp\;\n", @domains;
				$domains =~ s/<([c-zC-Z]).*?>//gs;
				$domains =~ s/<\/([c-zC-Z]).*?>//gs;
			}
			my $string; 
			if (length($domains)){
		  	$string = "$description<BR>Annotated domains:<BR>&nbsp;&nbsp;$domains";
			}
			else {
				$string = "$description"
			}
			WriteCache("$cache$val.ensembl.b",$string);
			return ($string || '');
		}
	}
	return $html;
}

sub isTb
{
	my ($_l) = @_;
	if($_l =~ /^AA[QZ]\d+$/ or $_l =~ /^CJA\d+$/ or $_l =~ /^EAN\d+$/)	{
		return 1;
	}
	return 0;
}
## GetCache opens the passed cached file name and returns the html from that file.
## called by GetInfo()

## updated as suggested by RH - 040721

sub GetCache
{
	my ($cache) = @_;
	my $cache = NewCache($cache);
	open(INPUT,"<$cache")  or return "";
	local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.

	my $return = <INPUT>;
	close INPUT;
	if($return =~ /www\.ebi\.uniprot\.org/)	{
		$return =~ s/www\.ebi\.uniprot\.org\/uniprot\-srv\/uniProtView\.do\?proteinAc\=/www\.uniprot\.org\/uniprot\//gsi;
	}
	return $return;
}

sub NewCache
{
	my ($_c) = @_;
	my $old = $_c;
	my ($v) = $_c =~ /.+\/(....)/;
	if(length($v) < 4)	{
		$v = "AAAA";
	}
	if($v =~ /\|/)	{
		($v) = $_c =~ /.+\/.+?\|(....)/;
	}
	$v =~ s/\|/_/g;
	if($old =~ /\.sp$/ or $old =~ /\.sp\.b$/)	{
		$v = "UPROT";
	}
	if($old =~ /\.genedb$/ or $old =~ /\.genedb\.b$/)	{
		$v = "GENEDB";
	}
	if($old =~ /\.toxo$/ or $old =~ /\.toxo\.b$/)	{
		$v = "TOXODB";
	}
	if(length($v) < 4)	{
		$v = "AAAA";
	}
	$_c =~ s/(.+\/)(.+?)/$1$v\/$2/;
	my ($d) = $_c =~ /(.+)\//;
	if(not -e $d)	{
		mkdir($d);
	}
	if(-e $old)	{
		my $c = "copy $old $_c";
		$c =~ s/\//\\/g;
		`$c`;
		unlink($old);
	}
	return $_c;
}

## WriteCache opens the passed cached file name and writes the passed html to that file.
## called by GetInfo()

## updated as suggested by RH - 040721

sub WriteCache
{
  	my ($val,$h) = @_;
	if(length($h) < 2)	{
		return;
	}
	return 0 unless (defined $val);  #actually not much reason to throw false when we aren't catching it anywhere...
  	$val = NewCache($val);
	open(INPUT,">$val") or return 0;
  	print INPUT $h;
  	close INPUT;
}


## DoSearch retieves and returns the html from the relevant database search site(ensmbl, ncbi etc..)
## by calling the appropriate Post() function. 
## called by homolog.pl, plist.pl and protein.pl.
## only the first two parameters are required: $l is the accession number and $protein is 0 or 1.
## if it is called from protein.pl, $protein is 1 and all the parameters are needed. this is so
## the from that holds the refresh button on that page contains the required hidden values.

sub DoSearch	{
	my ($l,$protein,$url,$label,$homolog,$uid_input,$refresh) = @_;
	$_ = $l;
	if(/\:reversed/)	{
		return "";
	}
	my $html = " ";
	if($protein == 1){
		$html = "<form method=\"POST\" action=\"/thegpm-cgi/protein.pl\">";
		$html .= "<input type=\"hidden\" name=\"path\" value=\"$url\" />";
		$html .= "<input type=\"hidden\" name=\"label\" value=\"$label\" />";
		$html .= "<input type=\"hidden\" name=\"homolog\" value=\"$homolog\" />";
		$html .= "<input type=\"hidden\" name=\"uid\" value=\"$uid_input\" />";
		$html .= "<input type=\"hidden\" name=\"refresh\" value=\"yes\" />";
		$html .= "<span class=\"small_label\">If the information below appears incomplete, press: </span><input type=\"submit\" class=\"but\" value=\"go\" /></form>\n";
	}
	if(/^HIT[0-9]+/)	{
		$html .= PostHit($l,$protein,$refresh);
		return $html;
	}	
	elsif(/^osa1[0-9]+/)	{
		$html .= PostTigrOsa($l,$protein,$refresh);
		return $html;
	}
	elsif(/IPI[0-9]+/ or /^IPI\:IPI[0-9]+/)	{
		$html .= PostIpi($l,$protein,$refresh);
		return $html;
	}	
	elsif(/gi\|/)	{
		$html .= PostNcbi($l,$protein,$refresh);
		return $html;
	}	
	elsif(/genedb\|/)	{
		$html .= PostGenedb($l,$protein,$refresh);
		return $html;
	}	
	elsif(/sp\|/)	{
		$html .= PostSProt($l,$protein,$refresh);
		return $html;
	}	
	elsif(/tr\|/)	{
		$html .= PostSProt($l,$protein,$refresh);
		return $html;
	}	
	elsif(/^At[1-9]g[0-9]+/)	{
		$html .= PostEnsemblMeta($l,$protein,$refresh);
		return $html;
	}	
	elsif(/^Bra[0-9]+\.\d/)	{
		$html .= PostEnsemblMeta($l,$protein,$refresh);
		return $html;
	}	
	elsif(/^tb[0-9]+/i)	{
		$html .= PostTigrTb($l,$protein,$refresh);
		return $html;
	}	
	elsif(/plasmoDB\|/i)	{
		$html .= PostPlasmoDB($l,$protein,$refresh);
		return $html;
	}
	elsif(/tgo\|/i)	{
		$html .= PostToxoDB($l,$protein,$refresh);
		return $html;
	}
	elsif(/^Y[A-Z]+?[0-9]+?[A-Z]/ or /^[RQ]\d{4}/)	{
		$html .= PostScd($l,$protein,$refresh);
#		$html .= PostMips($l,$protein,$refresh);
		return $html;
	}	
	elsif(/^AGAP/ or /^LOC_Os/ or /^AT[1-9C]G\d+/ or /[A-Z][0-9|A-Z]+?\.[0-9]/ or /^PPA\d+/ 
		or /^FBpp/ or /^DappuP/ or /^BRADI/ or /^SP[A-C][A-Z].*\-\d/ or /^CADAFUAP\d/ or isTb($_))	{
		$html .= PostEnsemblMeta($l,$protein,$refresh);
		return $html;
	}
	elsif(/DDB[0-9]+/ or /TA[0-9][0-9][0-9][0-9][0-9]/ or /SP[A-C][A-Z]/)	{
		$html .= PostDdb($l,$protein,$refresh);
		return $html;
	}	
	elsif(/ENS/ or /CG[0-9]+?-P[A-Z]/ or /GSTENP/ or /NEWSINFRUP/)	{
		$html .= PostEnsembl($l,$protein,$refresh);
		return $html;
	}
	return $html;
}


## returns the html from the relevant database search site and writes the new cache file
## called by DoSearch()

sub PostTigrAt
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
    	$l =~ s/\..+//;
	my $cache = get_cache_root() . "/cache/$l.kegg";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,60))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if(not($refresh eq "yes") and $ok == 0)	{
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
		return $html;
	}
	unlink($cache . ".b");
        my $url= "http://www.genome.jp/dbget-bin/www_bget?ath:$l";
	my @OutValue;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(60);
	my $req = GET($url);
	my $html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to KEGG.<BR>";
	}
	
	## updated to new regex style - RC 040721
	my ($def) = $html =~ /<body.+?\>(.*?)<\/body>/is;
	$def =~ s/(href=["|'])\//$1http\:\/\/www.genome.jp\//gs;	
	$def =~ s/(src=["|'])\//$1http\:\/\/www.genome.jp\//gs;	
	$def =~ s/<hr.*?$//si;
	$def =~ s/(spawn_help)/http:\/\/www.tigr.org\/tigr-scripts\/euk_manatee\/shared\/$1/;
	##
	
	my $value = qq(<script type="text/javascript">
<!-- Hide script
function btn(bobj,img) {
    bobj.src = "http://www.genome.jp/Fig/bget/button_" + img + ".gif";
}
function init(){
}
function Link_XtrctSeq2(form) {
    var dna_from;
    var dna_to;
    var dna_len;
    var plus_up   = Number(form.XtrctSeq_UP.value) ;
    var plus_down = Number(form.XtrctSeq_DOWN.value);
    var vector    = Number(form.VECTOR.value);
    var org       = form.ORG.value;
    var chr       = form.CHR.value;
    var kid       = form.KEGGID.value;
    var url;

    if (plus_up == 0 && plus_down == 0) {
	url = "http://www.genome.jp/dbget-bin/www_bget?-f+-n+n+" + kid;
    }
    else {
	if (vector == 1) {
	    dna_from  = Number(form.FROM.value)  - plus_up;
	    dna_to    = Number(form.TO.value) + plus_down;
	} else {
	    dna_from  = Number(form.FROM.value)  - plus_down;
	    dna_to    = Number(form.TO.value) + plus_up;
	}

	url = "/kegg-bin/cut_sequence_genes.pl?FROM=" + dna_from + "&TO=" + dna_to +"&VECTOR=" + vector + "&ORG=" + org;
	if (chr) url += "&CHR=" + chr;
    }
    //window.open( url, "_self" );
    location.href = url;
}
function go_taxonomy(form){
  form.submit();
}
// End script hiding --->
</script>

);
	$value .= "<hr><h3>KEGG protein report:</h3>" . $def;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostDdb
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $v = $l;
    	$l =~ s/\..+//;
	my $cache = get_cache_root() . "/cache/$l.ddb";
	$cache = NewCache($cache);
	my $ok = 0;
    	my $url= "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=dicty&name=$l&submit=Search";
	if($l =~ /TA[0-9][0-9][0-9][0-9][0-9]/)	{
    		$url= "http://www.genedb.org/genedb/Search?submit=Search+for&name=$l&organism=annulata&desc=yes&wildcard=yes";
		$cache = get_cache_root() . "/cache/$l.genedbta";
		$cache = NewCache($cache);
	}
	elsif($l =~ /^SP[A-C][A-Z]/)	{
    		$url= "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=pombe&name=$v&submit=Search";
		$cache = get_cache_root() . "/cache/$v.genedbsp";
		$cache = NewCache($cache);
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if(not($refresh eq "yes") and $ok == 0)	{
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
		return $html;
	}
	unlink($cache . ".b");
	my @OutValue;
	my $ref = ['formType'=>'navBar','name'=>$l,'organism' => 'dicty', 'submit' => 'Search'] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(60);
	my $req = GET($url);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to GeneDB.<BR>";
	}
	## updated to new regex style - RC 040721
	my ($def) = $html =~ /\<\!-- Main content part of page --\>(.*?)<\/body>/is;
	$def =~ s/(href=["|'])\//$1http\:\/\/www.genedb.org\//gs;	
	$def =~ s/(src=["|'])\//$1http\:\/\/www.genedb.org\//gs;	
	$def =~ s/<hr.*?$//si;
	##
	
	my $value = "<hr><h3>GeneDB protein report:</h3>" . "\n<script language=\"javascript\" src=\"http://www.genedb.org/zmenu.js\"></script>\n" . $def;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostTigrOsa
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $cache = get_cache_root() . "/cache/$l.tigrosa";
	$cache = NewCache($cache);
	my $ok = 0;
	open(INPUT,"<$cache")  or $ok = 1;
	if(not($refresh eq "yes") and $ok == 0)	{
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
		return $html;
	}
	unlink($cache . ".b");
    my $url= "http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi";
	my @OutValue;
	$l =~ s/osa1//;
	my $ref = ['db'=>'osa1','orf'=>$l] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));
	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(60);
	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to Tigr.<BR>";
	}
	if($html =~ /info is not currently available/)	{
		return "<BR><BR>Genomic information unavailable at Tigr.<BR>";
	}
	## updated to new regex style - RC 040721
	$html =~ s/.*?(\<table.*?\>)/$1/s;
	$html =~ s/\<hr.*?\>.*//si;
	$html =~ s/(href=["|'])\//$1http\:\/\/rice.plantbiology.msu.edu\//gsi;	
	$html =~ s/(src=["|'])\//$1http\:\/\/rice.plantbiology.msu.edu\//gsi;	
	$html =~ s/\<MAP.*\<\/MAP\>//si;
	my $value = "<hr><h3>MSU rice protein report:</h3>" . $html;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostTigrOsa5
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $cache = get_cache_root() . "/cache/$l.tigrosa5";
	$cache = NewCache($cache);
	my $ok = 0;
	open(INPUT,"<$cache")  or $ok = 1;
	if(not($refresh eq "yes") and $ok == 0)	{
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
		return $html;
	}
	unlink($cache . ".b");
    	my $url= "http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi";
	my @OutValue;
	$l =~ s/osa1//;
	my $ref = ['db'=>'osa1r5','orf'=>$l] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));
	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(60);
	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to Tigr.<BR>";
	}
	if($html =~ /info is not currently available/)	{
		return "<BR><BR>Genomic information unavailable at Tigr.<BR>";
	}
	## updated to new regex style - RC 040721
	$html =~ s/.*?(\<table.*?\>)/$1/s;
	$html =~ s/\<hr.*?\>.*//si;
	$html =~ s/(href=["|'])\//$1http\:\/\/rice.plantbiology.msu.edu\//gsi;	
	$html =~ s/(src=["|'])\//$1http\:\/\/rice.plantbiology.msu.edu\//gsi;	
	$html =~ s/\<MAP.*\<\/MAP\>//si;
	my $value = "<hr><h3>MSU rice protein report:</h3>" . $html;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostPlasmoDB
{
	my ($l,$protein,$refresh) = @_;
	return " ";
}

sub PostTigrTb
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $cache = get_cache_root() . "/cache/$l.tigrtb";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,60))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if(not($refresh eq "yes") and $ok == 0)	{
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
		return $html;
	}
	unlink($cache . ".b");
    my $url= "http://www.tigr.org/tigr-scripts/euk_manatee/shared/ORF_infopage.cgi";
	my @OutValue;
	my $ref = ['db'=>'tba1','orf'=>$l] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->timeout(10);
	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to Tigr.<BR>";
	}
	
	## updated to new regex style - RC 040721
	my ($def) = $html =~ /START OF MAIN DOCUMENT -->(.*?)<\/body>/is;
	$def =~ s/(href=["|'])\//$1http\:\/\/www.tigr.org\//gs;	
	$def =~ s/(src=["|'])\//$1http\:\/\/www.tigr.org\//gs;	
	$def =~ s/<hr.*?$//si;
	$def =~ s/(spawn_help)/http:\/\/www.tigr.org\/tigr-scripts\/euk_manatee\/shared\/$1/;
	
	my $value = "<hr><h3>Tigr protein report:</h3>" . $def;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostToxoDB
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	$l =~ s/tgo\|(.+?)\|/$1/;
	my $cache = get_cache_root() . "/cache/$l.toxo";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,60))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if(not($refresh eq "yes") and $ok == 0)	{
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
		return $html;
	}
	unlink($cache . ".b");
        my $url= "http://www.toxodb.org/toxo/showRecord.do";
	my @OutValue;

	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->timeout(10);
	my $req = GET("$url?name=GeneRecordClasses.GeneRecordClass&primary_key=$l&project_id=ToxoDB");
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to ToxoDB.<BR>";
	}
	## updated to new regex style - RC 040721
	my ($def) = $html =~ /genomic context --\>(.*)\<\/body\>/is;
	my ($name) = $html =~ /title\>.+?\((.+?)\).*?\<\/title/is;
	$name =~ s/\b+/ /gs;
	$def =~ s/\<\/body\>//gsi;
	$def =~ s/\<\/html\>//gsi;
	$def =~ s/(href\=["'])\//$1http\:\/\/www.toxodb.org\//gsi;
	$def =~ s/(src\=["'])\//$1http\:\/\/www.toxodb.org\//gsi;
	$def =~ s/(href\=["'])([a-gik-z])/$1http\:\/\/www.toxodb.org\/toxo\/$2/gsi;
	my ($js) = $html =~ /(\<script type\=\"text\/javascript\" src\=\'\/toxo\/js\/api\.js.+?'\>\<\/script\>)/si;
	$js =~ s/src\=\'/src=\'http:\/\/www.toxodb.org/s;
	$def = "<h3>$name</h3>" . $js . $def;
	my $value = "<hr><h3>ToxoDb protein report:</h3>" . $def;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}


## returns the html from the relevant database search site or cache if it available 
## and writes the new cache file
## called by DoSearch()
sub GetEnsemblArchive
{
	my ($_u) = @_;
	my ($l) = $_u =~ /(ENS[A-Z0-9]+)/;
	my ($old) = $_u =~ /\:\/\/(.+?)\//; 
	$old =~ tr/[A-Z]/[a-z]/;
	$_u =~ s/\:\/\/.+?\//\:\/\/$old\//;
	my $html;
	my $cache = get_cache_root() . "/cache/$l.ensembl";
	$cache = NewCache($cache);
	my $note;
	my $url= "$_u&show=snps&number=on";
    	my $local = "\"http:\/\/www.ensembl.org\/Homo_sapiens\/";
	my @OutValue;
	my $ref = ['peptide' => $l, 'show' => 'snps', 'number' => 'on'] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/5.0 (compatible; MSIE 7.0; Windows NT; DigExt)');
	$agent->timeout(20);
	my $req = GET($url);
#	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /\<body/si) or not($html =~ /Ensembl.*?\<\/th\>/si or $html =~ /The document has moved/))	{
		return "$url <br> $html";
	}
	##update to handle archived identifiers (eg: ENSP00000228652)
	my $retired = 0;
	my ($lower) = $_u =~ /(http\:\/\/.+?)\//;
	if($html =~ /\<b\>Retired\<\/b\>/si){
		($_) = $html =~ /\<a title\=\"View in archived protview\" href\=\"(.+?)\"/si;
		$ref = ['peptide' => $l, 'show' => 'snps', 'number' => 'on'] ;
		s/\?.+//;
		($lower) = $_ =~ /http\:\/\/(.+?)\//i;
		$lower =~ tr/[A-Z]/[a-z]/;
		s/http\:\/\/.+?\//http\:\/\/$lower\//;
		$req = $agent->post($_,$ref);
		$html = $req->content;
		$retired = 1;
		if(not($html =~ /\<body/si) or not($html =~ /Ensembl.*?\<\/th\>/si or $html =~ /The document has moved/))	{
			return "<BR><BR>Genomic information unavailable: cannot connect to ensembl. new url = \"$_\"<BR>";
			
		}
	}
	## updated to new regex style - RC 040721
	my ($def) = $html =~ /.*(<h3>.*? Protein Report<\/h3>.+)<\/body>/s;
	if(not length($def))	{
		($def) = $html =~ /.*(\<div id=\"page\"\>.*?)\<div class=\"sp\"\>/s;
	}
	$def =~ s/href=\"/href=\"$lower/isg;
	$def =~ s/src=\"/src=\"$lower/isg;
	$def =~ s/href=\"http:\/\/www.ensembl.orghttp:/href=\"http:/sg;
	$def =~ s/\<!-- begin footer --\>.*\<!-- end footer --\>//s;
	$def =~ s/\<form.*?\<\/form\>//s;
	$local .= "domainview\?domainentry\=";
	$def =~ s/\"domainview\?domainentry\=/$local/sg;
	if(length($def) > 2)	{
		open(INPUT,">$cache");
		print INPUT $def;
		close(INPUT);
	}
	return $note . $def;
}

sub PostEnsemblArchive
{
	my $html = "";
	my $note = "";
	my $acc;
	my @accs=();
	my ($l,$protein,$refresh) = @_;	
	my $cache = get_cache_root() . "/cache/$l.ensembl";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,60))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if($refresh ne "yes" and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			## check that the label is not one of the ones that ensembl has replaced
			if(not $html =~ /Identifier Removed from Database/si){
				return $note . $html;
			} 
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    my $url= "http://www.ensembl.org/Homo_sapiens/protview";
    my $local = "\"http:\/\/www.ensembl.org\/Homo_sapiens\/";
    ##added 20050215
    my ($label) = $l =~ /(ENS\w*[PG])\d+/; 
    if($label){
		$url = "http://www.ensembl.org/$g_ens{$label}/protview";
		$local = "\"http:\/\/www.ensembl.org\/$g_ens{$label}\/";
	}
	elsif($l =~ /[A-Z][0-9|A-Z]+?\.[0-9]/)	{
		$url = "http://www.ensembl.org/Caenorhabditis_elegans/protview";
		$local = "\"http:\/\/www.ensembl.org\/Caenorhabditis_elegans\/";
	}
	elsif($l =~ /CG[0-9]+?-P[A-Z]/)	{
		$url = "http://www.ensembl.org/Drosophila_melanogaster/protview";
		$local = "\"http:\/\/www.ensembl.org\/Drosophila_melanogaster\/";
	}
	elsif($l =~ /^FBpp/)	{
		$url = "http://www.ensembl.org/Drosophila_melanogaster/protview";
		$local = "\"http:\/\/www.ensembl.org\/Drosophila_melanogaster\/";
	}
	elsif($l =~ /^DappuP/)	{
		$url = "http://www.ensembl.org/Daphnia_pulex/protview";
		$local = "\"http:\/\/www.ensembl.org\/Daphnia_pulex\/";
	}
	elsif($l =~ /^AGAP/)	{
		$url = "http://www.ensembl.org/Anopheles_gambiae/protview";
		$local = "\"http:\/\/www.ensembl.org\/Anopheles_gambiae\/";
	}
	elsif($l =~ /NEWSINFRUP/)	{
		$label = "NEWSINFRUP";
		$url = "http://www.ensembl.org/$g_ens{$label}/protview";
		$local = "\"http:\/\/www.ensembl.org\/$g_ens{$label}\/";
	}
	elsif($l =~ /GSTENP/)	{
		$label = "GSTENP";
		$url = "http://www.ensembl.org/$g_ens{$label}/protview";
		$local = "\"http:\/\/www.ensembl.org\/$g_ens{$label}\/";
	}
	my @OutValue;
	my $ref = ['peptide' => $l, 'show' => 'snps', 'number' => 'on'] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/5.0 (compatible; MSIE 7.0; Windows NT; DigExt)');
	$agent->timeout(20);
	my $req = GET("$url?peptide=$l&show=snps&number=on");
#	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /\<body/si) or not($html =~ /Ensembl.*?\<\/th\>/si or $html =~ /The document has moved/))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to ensembl.<BR>";
	}
	##update to handle archived identifiers (eg: ENSP00000228652)
	if($html =~ /The document has moved/){
		($_) = $html =~ /href=\"(.+?)\"/s;
		$req = GET("http://www.ensembl.org/$_");
		$html = $agent->request($req)->as_string();
		if($html =~ /This ID has been removed from Ensembl/s)	{
			($_) = $html =~ /This ID has been removed from Ensembl.+?\<a .+?\>(.+?)\<\/a\>/s;
			$ref = ['peptide' => $_, 'show' => 'snps', 'number' => 'on'] ;
			$req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
			$html = $agent->request($req)->as_string();
		}
		if(not($html =~ /\<body/si) or not($html =~ /Ensembl.*?\<\/th\>/si or $html =~ /The document has moved/))	{
			return "<BR><BR>Genomic information unavailable: cannot connect to ensembl.<BR>";
		}
		s/.+?\=//;
		if($l ne $_)	{
			$note = "<span class=\"alt\">NOTE 1: Ensembl has replaced the accession number: '$l' with '$_'.</span>";
		}
	} 
	my $retired = 0;
	my $lower;
	if($html =~ /\<b\>Retired\<\/b\>/si){
		($_) = $html =~ /\<a title\=\"View in archived protview\" href\=\"(.+?)\"/si;
		$ref = ['peptide' => $l, 'show' => 'snps', 'number' => 'on'] ;
		s/\?.+//;
		($lower) = $_ =~ /http\:\/\/(.+?)\//i;
		$lower =~ tr/[A-Z]/[a-z]/;
		s/http\:\/\/.+?\//http\:\/\/$lower\//;
		$req = $agent->post($_,$ref);
		$html = $req->content;
		$retired = 1;
		if(not($html =~ /\<body/si) or not($html =~ /Ensembl.*?\<\/th\>/si or $html =~ /The document has moved/))	{
			return "<BR><BR>Genomic information unavailable: cannot connect to ensembl. new url = \"$_\"<BR>";
			
		}
	}
	## updated to new regex style - RC 040721
	my ($def) = $html =~ /.*(<h3>.*? Protein Report<\/h3>.+)<\/body>/s;
	if(not length($def))	{
		($def) = $html =~ /.*(\<div id=\"page\"\>.*?)\<div class=\"sp\"\>/s;
	}
	if($retired)	{
		$def =~ s/href=\"/href=\"http:\/\/$lower/isg;
		$def =~ s/src=\"/src=\"http:\/\/$lower/isg;
	}
	else	{
		$def =~ s/href=\"/href=\"http:\/\/www.ensembl.org/isg;
		$def =~ s/src=\"/src=\"http:\/\/www.ensembl.org/isg;
	}
	$def =~ s/href=\"http:\/\/www.ensembl.orghttp:/href=\"http:/sg;
	$def =~ s/\<!-- begin footer --\>.*\<!-- end footer --\>//s;
	$def =~ s/\<form.*?\<\/form\>//s;
	$local .= "domainview\?domainentry\=";
	$def =~ s/\"domainview\?domainentry\=/$local/sg;
	if(length($def) > 2)	{
		open(INPUT,">$cache");
		print INPUT $def;
		close(INPUT);
	}
	return $note . $def;
}

sub PostEnsembl
{
	my $html = "";
	my $note = "";
	my $acc;
	my @accs=();
	my ($l,$protein,$refresh) = @_;	
	my $cache = get_cache_root() . "/cache/$l.ensembl";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,5))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if(not $ok){
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
	}
	if(not $html =~ /Description/si)	{
		$protein = 0;
		$ok = 1;
	}
	if($html =~ /Unable to produce objects \- panic\!/si)	{
	}
	elsif($protein == 1){
		if($refresh ne "yes" and $ok == 0)	{
			## check that the label is not one of the ones that ensembl has replaced
			if(not $html =~ /Identifier Removed from Database/si){
				return $note . $html;
			} 
		}
	}
	elsif($ok == 0)	{
		return $html;
	}
	unlink($cache . ".b");
    my $url= "http://www.ensembl.org/Homo_sapiens/protview";
    my $local = "\"http:\/\/www.ensembl.org\/Homo_sapiens\/";
    ##added 20050215
    my ($label) = $l =~ /(ENS\w*[PG])\d+/; 
    my $variation = "exons=yes;";
    if($label){
		$url = "http://www.ensembl.org/$g_ens{$label}/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/$g_ens{$label}\/";
		if($label eq "ENSP" or $label eq "ENSMUSP" or $label eq "ENSRNOP")	{
			$variation = "exons=yes;variation=yes;";
		}
	}
	elsif($l =~ /[A-Z][0-9|A-Z]+?\.[0-9]/)	{
		$url = "http://www.ensembl.org/Caenorhabditis_elegans/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/Caenorhabditis_elegans\/";
		$variation = "exons=yes;variation=yes;";
	}
	elsif($l =~ /CG[0-9]+?-P[A-Z]/)	{
		$url = "http://www.ensembl.org/Drosophila_melanogaster/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/Drosophila_melanogaster\/";
	}
	elsif($l =~ /^DappuP/)	{
		$url = "http://www.ensembl.org/Daphnia_pulex/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/Daphnia_pulex\/";
		$label = "DappuP";
		$variation = "exons=yes;variation=yes;";
	}
	elsif($l =~ /^FBpp/)	{
		$url = "http://www.ensembl.org/Drosophila_melanogaster/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/Drosophila_melanogaster\/";
		$label = "FBpp";
		$variation = "exons=yes;variation=yes;";
	}
	elsif($l =~ /^AGAP/)	{
		$url = "http://www.ensembl.org/Anopheles_gambiae/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/Anopheles_gambiae\/";
		$label = "AGAP";
	}
	elsif($l =~ /NEWSINFRUP/)	{
		$label = "NEWSINFRUP";
		$url = "http://www.ensembl.org/$g_ens{$label}/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/$g_ens{$label}\/";
	}
	elsif($l =~ /GSTENP/)	{
		$label = "GSTENP";
		$url = "http://www.ensembl.org/$g_ens{$label}/Transcript/ProteinSummary";
		$local = "\"http:\/\/www.ensembl.org\/$g_ens{$label}\/";
	}
####
######## New Ensembl Interface ########
####
	my $ensembl = 'http://www.ensembl.org/';  
	my $species = $g_ens{$label};
	if(not $species)	{
		$species = 'Caenorhabditis_elegans';
	}
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	$agent->agent('Mozilla/4.0 (compatible; MSIE 6.0; Windows NT; DigExt)');
	my $timeout = 20;
	$agent->timeout($timeout);
	$url .= "?p=$l";
	my $req = new HTTP::Request GET => $url;
	my $res = $agent->request($req);
	my $html;
	if ($res->is_success) {
		$html = $res->content;
	} 
	else	{
		$timeout = 40;
		$agent->timeout($timeout);
		$req = new HTTP::Request GET => $url;
		$res = $agent->request($req);
		$html = $res->content;
		if ($res->is_success) {
			$html = $res->content;
		}
		else	{ 
			return "<!-- Ensembl not available: url = $url, html = $html-->";
		}
	}
	if($html =~ /This transcript is not in the current gene set/) {
		my ($u) = $html =~ /Latest version.+?\<a href\=\"(.+?)\"/si;
		my ($v) = $u =~ /\/\/(.+?)\//;
		$v =~ tr/[A-Z]/[a-z]/;
		$url =~ s/www\.ensembl\.org/$v/;
		$ensembl = "http://$v/";
		$timeout = 40;
 		$agent->timeout($timeout);
		$req = new HTTP::Request GET => $url;
		$res = $agent->request($req);
		$html = $res->content;
		if ($res->is_success) {
			$html = $res->content;
		} 
		else	{
			return "<!-- Ensembl not available: url = $url, html = $html -->";
		}
		##GetEnsemblArchive($u);
	}
	my ($link) = $html =~ /\<a href\=\".+\?(.+?)\"\>$l\<\/a\>/si;
	my $desc = $html;
	$desc =~ s/.+?\<div class\=\"content\"\>//si;
	if($desc =~ /\"logo_holder\"/si)	{
		$desc =~ s/.+?\<div class\=\"content\"\>//si;
	}
	$desc =~ s/&lt;/\</gsi;
	$desc =~ s/&gt;/\>/gsi;
	$desc =~ s/&quot;/\"/gsi;
	$desc =~ s/\<\/div.+/<\/div>/si; 
	$desc =~ s/src\=\"\//src\=\"$ensembl/gi;
	$desc =~ s/href\=\"\//href\=\"$ensembl/gi;
	$desc =~ s/href *\= *\//href \= $ensembl\//gi;
	$desc =~ s/( href *\=)/ target="_ENSEMBL"$1/gi;
	$desc =~ s/\<\/table\>.+/\<\/table\>\<\/div\>/si;
	my %modules = (	"DomainSpreadsheet"	=> "" ,
		"TranslationImage"	=> "" ,
		"ProteinVariations"	=> "" ,
		"ProteinSeq"		=> "" ,
		"SimilarityMatches"	=> "");
	$desc =~ s/\<p\>/<p>Description: /;
	if($l =~ /^FBpp/)	{
		$desc = "<p>Description: no description available</p>" . $desc;
	}
	my $width = 700;
	my $value = "TranslationImage";
	$html = qq(<!-ENSEMBL60-->
	<table width="600">
	<tr>
		<td><b>ENSEMBL Protein report: $l</b></td>
	</tr>
	<tr>
		<td><BR>$desc</td>
	</tr>
	);
	$html .= "<table width=\"$width\"><tr>\t<td>";
	my $module = getModule($link, $ensembl, $value, $species, $timeout,$agent);
	my ($png) = $module =~ /\<a ([^\>]+?)\>Export as PNG\<\/a\>/smi;
	($png) = $png =~ /href\=\"(.+?)\"/;
	if(not $png)	{
		($png) = $local =~ /^\"(.+)/;
		$png .= "Component/Transcript/Web/TranslationImage?db=core;t=$l;export=png";
	}
	my $oens = getPng($png,$l,$timeout,$agent);
	$module =~ s/src\=\".+?\"/src\=\"$oens\"/smi;
	$module =~ s/(<img .+?\>).+/$1\n\<\/div\>\<\/div\>/smi;
	$html .= $module;
	$html .= "\t</td>\n</tr>\n</table>\n";
	$value = "ProteinSeq";
	$html .= "<table width=\"$width\"><tr>\t<td valign=\"top\"><b>Protein Sequence:</b></td><td valign=\"top\">";
	my $temp = getModule("$link;$variation", $ensembl, $value, $species, $timeout,$agent);
	$temp =~ s/\<pre\>/<pre>\n/m;
	$html .= $temp;
	$html .= "\t<hr></td>\n</tr>\n</table>\n";
	$html .= "<!-- ProteinSeq url = $link;exons=yes;variation=yes; -->\n";
	$value = "DomainSpreadsheet";
	$html .= "<table width=\"$width\"><tr>\t<td>";
#	$html .= "<hr>$link, $ensembl, $value, $species, $timeout<hr>";
	$html .= getModule($link, $ensembl, $value, $species, $timeout,$agent);
	$html .= "\t</td>\n</tr>\n</table>\n";
	$value = "ProteinVariations";
	$html .= "<table width=\"$width\"><tr>\t<td>";
	$html .= getModule($link, $ensembl, $value, $species, $timeout,$agent);
	$html .= "\t</td>\n</tr>\n</table>\n";
	$value = "SimilarityMatches";
	$html .= "<table width=\"$width\"><tr>\t<td>";
	my $sim = getModule($link, $ensembl, $value, $species, $timeout,$agent,$agent);
	if($l =~ /^FBpp/ and $sim =~ /NP\_\d+/)	{
		my $d = $sim;
		$d =~ s/.+?NP\_\d+\.\d//s;
		$d =~ s/.+?align\<\/a\>] \<br \/\>//s;
		$d =~ s/\[.+//s;
		$html =~ s/Description: no description available/Description: $d/s;
	}
	$html .= $sim;
	$html .= "\t</td>\n</tr>\n";
	$html .= "</table>\n";
	if(length($html) > 2)	{
		open(INPUT,">$cache");
		print INPUT $html;
		close(INPUT);
	}
	return $html;
}

sub PostEnsemblMeta
{
	my $html = "";
	my $note = "";
	my $acc;
	my @accs=();
	my ($l,$protein,$refresh) = @_;	
	my $cache = get_cache_root() . "/cache/$l.ensembl";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,60))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if(not $ok){
		local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
		$html = <INPUT>;
		close(INPUT);
	}
	if(not $html =~ /Description/si)	{
		$protein = 0;
		$ok = 1;
	}
	if($html =~ /Unable to produce objects \- panic\!/si)	{
	}
	elsif($protein == 1){
		if($refresh ne "yes" and $ok == 0)	{
			## check that the label is not one of the ones that ensembl has replaced
			if(not $html =~ /Identifier Removed from Database/si){
				return $note . $html;
			} 
		}
	}
	elsif($ok == 0)	{
		return $html;
	}
	unlink($cache . ".b");
    	my $url= "http://www.ensembl.org/Homo_sapiens/protview";
    	my $local = "\"http:\/\/www.ensembl.org\/Homo_sapiens\/";
    ##added 20050215
	my $kingdom = 'metazoa';
    	my ($label) = $l =~ /([A-Z\_]+)\d+/i; 
	my $species = $g_meta{$label};
	my $alt = $l;
	my $image;
     	if($label =~ /LOC\_Os/){
		$url = "http://plants.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;t=";
		$local = "\"http:\/\/plants.ensembl.org\/$g_ens{$label}\/";
		$kingdom = 'plants';
	     }
     	elsif($label =~ /^AT/i){
		$label =~ tr/[a-z]/[A-Z]/;
		$url = "http://plants.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;p=";
		$local = "\"http:\/\/plants.ensembl.org\/$g_meta{$label}\/";
		if($l !~ /\.\d/)	{
			$l .= ".1";
			$alt = $l;
		}
		$alt =~ s/\-P/\-TAIR/;
		$species = 'Arabidopsis_thaliana';
		$kingdom = 'plants';
     	}
      	elsif($label =~ /^Bra/i){
		$url = "http://plants.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;p=";
		$local = "\"http:\/\/plants.ensembl.org\/$g_meta{$label}\/";
		if($l !~ /\.\d/)	{
			$l .= ".1";
			$alt = $l;
		}
		$species = $g_meta{$label};
		$kingdom = 'plants';
     	}
      	elsif($label =~ /^AA[QZ]$/i or $label =~ /^EAN$/ or $label =~ /^CAJ$/){
		$url = "http://protists.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;p=";
		$local = "\"http:\/\/protists.ensembl.org\/$g_meta{$label}\/";
		$species = $g_meta{$label};
		$alt = $l;
		$kingdom = 'protists';
     	}
    	elsif($label =~ /BRADI/){
		$url = "http://plants.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;t=";
		$local = "\"http:\/\/plants.ensembl.org\/$g_ens{$label}\/";
		$alt =~ s/\-P/\-TAIR/;
		$kingdom = 'plants';
     	}
	elsif($l =~ /SP[A-Z]+.+\-\d/)	{
		$url = "http://fungi.ensembl.org/Schizosaccharomyces_pombe/Transcript/ProteinSummary?db=core;t=";
		$local = "\"http:\/\/www.ensembl.org\/Schizosaccharomyces_pombe\/";
		$species = "Schizosaccharomyces_pombe";
		$kingdom = "fungi";
	}
	elsif($l =~ /CADAFUAP\d/)	{
		$url = "http://fungi.ensembl.org/Aspergillus_fumigatus/Transcript/ProteinSummary?db=core;t=";
		$local = "\"http:\/\/fungi.ensembl.org\/Aspergillus_fumigatus\/";
		$species = "Aspergillus_fumigatus";
		$kingdom = "fungi";
	}
	elsif($l =~ /[A-Z][0-9|A-Z]+?\.[0-9]/)	{
		$url = "http://metazoa.ensembl.org/Caenorhabditis_elegans/Transcript/ProteinSummary?db=core;t=";
		$local = "\"http:\/\/www.ensembl.org\/Caenorhabditis_elegans\/";
		$species = "Caenorhabditis_elegans";
		$kingdom = 'metazoa';
	}
    	elsif($label =~ /DappuP/){
		$url = "http://metazoa.ensembl.org/$g_meta{$label}/protview?peptide=";
		$local = "\"http:\/\/metazoa.ensembl.org\/$g_ens{$label}\/";
 		$kingdom = 'metazoa';
   	}
    	elsif($label =~ /FBpp/){
		$url = "http://metazoa.ensembl.org/$g_meta{$label}/protview?peptide=";
		$local = "\"http:\/\/metazoa.ensembl.org\/$g_ens{$label}\/";
 		$kingdom = 'metazoa';
   	}
     	elsif($label =~ /PPA/){
		$url = "http://metazoa.ensembl.org/Pristionchus_pacificus/Transcript/ProteinSummary?db=core;t=";
		$local = "\"http:\/\/metazoa.ensembl.org\/Pristionchus_pacificus\/";
 		$species = "Pristionchus_pacificus";
		$kingdom = 'metazoa';
   	}
   	else{
		$url = "http://metazoa.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;t=";
		$local = "\"http:\/\/metazoa.ensembl.org\/$g_ens{$label}\/";
 		$kingdom = 'metazoa';
   	}

####
######## New Ensembl Interface ########
####
	my $ensembl = "http://$kingdom.ensembl.org/";  
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	$agent->agent('Mozilla/4.0 (compatible; MSIE 6.0; Windows NT; DigExt)');
	my $timeout = 40;
	$agent->timeout($timeout);
	$url .= "$alt";
	my $req = new HTTP::Request GET => $url;
	my $res = $agent->request($req);
	my $html;
	if ($res->is_success) {
		$html = $res->content;
	} 
	else	{
		$timeout = 60;
		$agent->timeout($timeout);
		$req = new HTTP::Request GET => $url;
		$res = $agent->request($req);
		if ($res->is_success) {
			$html = $res->content;
		}
		else	{ 
			return "<!-- Ensembl not available: url = $url, html = $html-->";
		}
	}
	if($html =~ /This transcript is not in the current gene set/) {
		my ($u) = $html =~ /Latest version.+?\<a href\=\"(.+?)\"/si;
		my ($v) = $u =~ /\/\/(.+?)\//;
		$v =~ tr/[A-Z]/[a-z]/;
		$url =~ s/$kingdom\.ensembl\.org/$v/;
		$ensembl = "http://$v/";
		$timeout = 60;
 		$agent->timeout($timeout);
		$req = new HTTP::Request GET => $url;
		$res = $agent->request($req);
		if ($res->is_success) {
			$html = $res->content;
		} 
		else	{
			return "<!-- Ensembl not available -->";
		}
		##GetEnsemblArchive($u);
	}
	my ($link) = $html =~ /\<a href\=\".+\?(.+?)\"\>$l\<\/a\>/si;
	if(!$link)	{
		($link) = $html =~ /\<a href\=\".+\?(.+?)\"\>$l\-P\<\/a\>/si;
	}
	my $desc = $html;
	$desc =~ s/.+?\<div class\=\"content\"\>//si;
	if($desc =~ /\"logo_holder\"/si)	{
		$desc =~ s/.+?\<div class\=\"content\"\>//si;
	}
	$desc =~ s/&lt;/\</gsi;
	$desc =~ s/&gt;/\>/gsi;
	$desc =~ s/&quot;/\"/gsi;
	$desc =~ s/\<\/div.+/<\/div>/si; 
	$desc =~ s/src\=\"\//src\=\"$ensembl/gi;
	$desc =~ s/href\=\"\//href\=\"$ensembl/gi;
	$desc =~ s/href *\= *\//href \= $ensembl\//gi;
	$desc =~ s/( href *\=)/ target="_ENSEMBL"$1/gi;

	my %modules = (	"DomainSpreadsheet"	=> "" ,
		"TranslationImage"	=> "" ,
		"ProteinVariations"	=> "" ,
		"ProteinSeq"		=> "" ,
		"SimilarityMatches"	=> "");
	$desc =~ s/\<p\>/<p>Description: /;
	if($l =~ /^FBpp/)	{
		$desc = "<p>Description: no description available</p>" . $desc;
	}
	my $width = 700;
	my $value = "TranslationImage";
	$html = qq(<!-ENSEMBL60-->
	<table width="600">
	<tr>
		<td><b>ENSEMBL Protein report: $l</b></td> 
	</tr>
	<tr>
		<td><BR>$desc</td>
	</tr>
	);
	$html .= "<table width=\"$width\"><tr>\t<td>";
	my $module = getModule($link, $ensembl, $value, $species, $timeout,$agent);
	my ($png) = $module =~ /\<a ([^\>]+?)\>Export as PNG\<\/a\>/smi;
	($png) = $png =~ /href\=\"(.+?)\"/;
	if(not $png)	{
		($png) = $local =~ /^\"(.+)/;
		$png .= "Component/Transcript/Web/TranslationImage?db=core;t=$l;export=png";
	}
	my $oens = getPng($png,$l,$timeout,$agent);
	if($oens !~ /\<\!/)	{
		$module =~ s/src\=\".+?\"/src\=\"$oens\"/smi;
		$module =~ s/(<img .+?\>).+/$1\n\<\/div\>\<\/div\>/smi;

	}
	else	{
		$module .= "\n\n $oens \n\n";
	}
	$html .= $module;
	$html .= "\t</td>\n</tr>\n</table>\n";
	$value = "ProteinSeq";
	$html .= "<table width=\"$width\"><tr>\t<td valign=\"top\"><b>Protein Sequence:</b></td><td valign=\"top\">";
	my $temp;
	if(0 and $label !~ /^PPA/)	{
		$temp = getModule("$link;exons=yes;variation=yes;", $ensembl, $value, $species, $timeout,$agent);
	}
	else	{
		$temp = getModule("$link;exons=yes;", $ensembl, $value, $species, $timeout,$agent);
	}
	$temp =~ s/\<pre\>/<pre>\n/m;
	$html .= $temp;
	$html .= "\t<hr></td>\n</tr>\n</table>\n";
	$value = "DomainSpreadsheet";
	$html .= "<table width=\"$width\"><tr>\t<td>";
#	$html .= "<hr>$link, $ensembl, $value, $species, $timeout<hr>";
	$html .= getModule($link, $ensembl, $value, $species, $timeout,$agent);
	$html .= "\t</td>\n</tr>\n</table>\n";
	$value = "ProteinVariations";
	$html .= "<table width=\"$width\"><tr>\t<td>";
	$html .= getModule($link, $ensembl, $value, $species, $timeout,$agent);
	$html .= "\t</td>\n</tr>\n</table>\n";
	$value = "SimilarityMatches";
	$html .= "<table width=\"$width\"><tr>\t<td>";
	my $sim = getModule($link, $ensembl, $value, $species, $timeout,$agent);
	if($l =~ /^FBpp/ and $sim =~ /NP\_\d+/)	{
		my $d = $sim;
		$d =~ s/.+?NP\_\d+\.\d//s;
		$d =~ s/.+?align\<\/a\>] \<br \/\>//s;
		$d =~ s/\[.+//s;
		$html =~ s/Description: no description available/Description: $d/s;
	}
	$html .= $sim;
	$html .= "\t</td>\n</tr>\n";
	$html .= "</table>\n";
	if(length($html) > 2)	{
		open(INPUT,">$cache");
		print INPUT $html;
		close(INPUT);
	}
	return $html;
}

sub getModule
{
	my ($_link,$_ens,$_id,$_sp,$_to,$_ua) = @_;
	my $route = "/$_sp/Component/Transcript/Web/";
	my $url = "$_ens$route$_id?$_link";
	my $agent;
	if($_ua)	{
		$agent = $_ua;
	}
	else	{
		$agent = LWP::UserAgent->new();
		$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));
		$agent->agent('Mozilla/4.0 (compatible; MSIE 6.0; Windows NT; DigExt)');
	}
	if($_to)	{
		$agent->timeout($_to);
	}
	else	{
		$agent->timeout(10);
	}
	my $req = new HTTP::Request GET => $url;
	my $res = $agent->request($req);
 	my $html;
	if ($res->is_success)	{
  		$html = "\n\n<!-- Ensembl OK, url=$url -->\n\n" . $res->content;
	} 
	else	{
 	 	return "<!-- Ensembl not available, url=$url -->";
	}
	$html =~ s/src\=\"\//src\=\"$_ens/gi;
	$html =~ s/href\=\"\//href\=\"$_ens/gi;
	$html =~ s/href *\= *\//href \= $_ens\//gi;
	$html =~ s/( href *\=)/ target="_ENSEMBL"$1/gi;
	$html =~ s/<h2>/<b>/g;
	$html =~ s/<\/h2>/<\/b>/g;
	$html =~ s/width\:100\%;/width\:75\%;/gi;
	return $html;
}

sub getPng
{
	my ($_u,$_ens,$_to,$_ua) = @_;
	my $ag;
	if($_ua)	{
		$ag = $_ua;
	}
	else	{
		$ag = LWP::UserAgent->new();
		$ag->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));
		$ag->agent('Mozilla/4.0 (compatible; MSIE 6.0; Windows NT; DigExt)');
	}
	$_to = 60;
	if($_to)	{
		$ag->timeout($_to);
	}
	else	{
		$ag->timeout(30);
	}
	my $rq = new HTTP::Request GET => $_u;
	my $rs = $ag->request($rq);
 	if ($rs->is_success)	{
		my $p = "/llama-magic/archive/pics/$_ens.png";
		open(OUT,">..$p");
		binmode OUT;
  		print OUT $rs->content;
		close(OUT);
		return $p;
 	} 
	else	{
		my $html = $rs->content;
 	 	return "<!-- Ensembl image not available: url = $_u, timeout = $_to, html = $html -->";
	}
}

####
######## End new Ensembl Interface ########
####

## returns the html from the relevant database search site or cache if it available 
## and writes the new cache file
## called by DoSearch()

sub PostNcbi
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $gi = $l;
	my $nuc = 0;
	$gi =~ s/gi\|//s;
	$gi =~ s/\|.*//s;
	my $cache = get_cache_root() . "/cache/$gi.ncbi";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,14))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			return $html;
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    	my $url= "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&db=protein&val=$gi&dopt=genpept&sendto=on&log\$=seqview";
	my @OutValue;
	if($nuc == 1)	{
                $url= "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&db=nucleotide&val=$gi&sendto=on&log\$=seqview";
        }
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(10);
	my $req = GET($url);
	$html = $agent->request($req)->as_string();
	$html =~ s/.+?LOCUS/\<pre\>\nLOCUS/si;
	$html = "$html</pre>\n";
	if($html =~ /GenPept does not exist for gi/si or $html =~ /does not have presentation/si or $html =~ /GenPept is not available/si or $html =~ /Requested format \'GenPept\' is not applicable/)	{
		$url= "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&db=nucleotide&val=$gi&sendto=on&log\$=seqview";
 		$req = GET($url);
		$html = $agent->request($req)->as_string();
	}
	if(not($html =~ /locus/si) or ($html =~ /WWW\sError/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to ncbi.<BR>";
	}
	## updated to new regex style - RC 040721
	$html =~ s/.*(<pre.*?>)/$1\n/si;
	$html =~ s/(<\/pre>).*/$1\n/si;
	$html =~ s/(href="?)\//$1http\:\/\/www.ncbi.nlm.nih.gov\//gs;
	my $value = "<hr><h3>NCBI protein report:</h3>" . $html;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostGenedb
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $gi = $l;
	my $nuc = 0;
	$gi =~ s/genedb\|//s;
	$gi =~ s/\|.*//s;
	my $cache = get_cache_root() . "/cache/$gi.genedb";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,1))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			return $html;
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    	my $url= "http://www.genedb.org/gene/$gi";
	my @OutValue;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(10);
	my $req = GET($url);
	$html = $agent->request($req)->as_string();
	$html =~ s/.+?(\<h2.*?\>[\&nbsp\; ]*?General Information)/$1/si;
	if(not($html =~ /General Information/si))	{
		return "<BR><BR>Genomic information unavailable from GeneDB.<BR>";
	}
	$html =~ s/(href="?)\//$1http\:\/\/www.genedb.org\//gs;
	## updated to new regex style - RC 040721
	my ($out) = $html =~ /(.+?\<\/table\>)/si;
	my $go = "";
	my $ppd = "";
	if($html =~ /\<h2.*?\>[\&nbsp\; ]*?Gene Ontology/si)	{
		$html =~ s/.+?(\<h2\>[\&nbsp\; ]*?Gene Ontology)/$1/si;
		($go) = $html =~ /(.+?\<\/table\>)/si;
	}
	if($html =~ /\<h2.*?\>[\&nbsp\; ]*?Predicted Peptide Data/si)	{
		$html =~ s/.+?(\<h2\>[\&nbsp\; ]*?Predicted Peptide Data)/$1/si;
		($ppd) = $html =~ /(.+?\<\/table\>)/si;
	}
	my $value = "<hr><h3>GeneDB protein report:</h3>" . "$out\n$go\n$ppd\n<BR><BR>\n";
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostSProt
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $gi = $l;
	if($gi =~ /^sp\|/)	{
		$gi =~ s/sp\|//s;
	}
	if($gi =~ /^tr\|/)	{
		$gi =~ s/tr\|//s;
	}
	$gi =~ s/\|.*//s;
	my $cache = get_cache_root() . "/cache/$gi.sp";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,60))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			return $html;
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    	my $url= "http://www.uniprot.org/uniprot/$gi.txt";
	my @OutValue;
	my $ref = ['id' => $gi] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(10);
	my $req = GET($url);
	$html = $agent->request($req)->as_string();
	if($html =~ /Location\: \/uniprot\//)	{
		($url) = $html =~ /Location: (\/uniprot\/.+?format\=txt)/;
		$url = "http://www.uniprot.org" . $url;
		$req = GET($url);
		$html = $agent->request($req)->as_string();
		my ($up) = $html =~ /\<(\W+)\.rdf/si;
		if($up)	{
			$url = "http://www.uniprot.org/uniprot/$up.txt";
			$req = GET($url);
			$html = $agent->request($req)->as_string();
		}
	}
	else	{
		my ($up) = $html =~ /(\w+)\.rdf/s;
		if($up)	{
			$url = "http://www.uniprot.org/uniprot/$up.txt";
			$req = GET($url);
			$html = $agent->request($req)->as_string();
		}
	}
	if(not $html =~ /ID /s)	{
		return "<BR><BR>Genomic information unavailable: cannot connect to UniProt at $url<BR>";
	}
	## updated to new regex style - RC 040721
	$html =~ s/.+?ID /ID /si;
	my $value = "<hr><pre>$html</pre>";
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}


sub PostIpi
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	if($l =~ /^IPI\:IPI[0-9]+\./)	{
		$l =~ s/^IPI\://;
	}
	$l =~ s/\.[0-9]*//;
	my $cache = get_cache_root() . "/cache/$l.ipi";
	$cache = NewCache($cache);
	my $ok = 0;
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			return $html;
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    my $url= "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-e+[IPI-acc:$l]+-vn+2";
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(10);
	my $req = GET($url);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<pre>.*<\/pre>/si))	{
		return "<BR><BR>Protein information unavailable: cannot connect to IPI.<BR>";
	}
	## updated to new regex style - RC 040721
	$html =~ s/.*(<pre>)/$1\n/si;
	$html =~ s/(<\/pre>).*/$1\n/si;
	$html =~ s/href=\"(wgetz)/href=\"http:\/\/srs.ebi.ac.uk\/srsbin\/cgi-bin\/$1/gsi;	
	my $value = "<hr><h3>IPI protein report:</h3>" . $html;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}


## returns the html from the relevant database search site or cache if it available 
## and writes the new cache file
## called by DoSearch()

sub PostWorm
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $cache = get_cache_root() . "/cache/$l.worm";
	$cache = NewCache($cache);
	my $ok = 0;
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			return $html;
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    my $url= "http://legacy.wormbase.org/db/gene/gene";
	my @OutValue;
	my $ref = ['name' => $l] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(10);
	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to wormbase.<BR>";
	}
	## updated to new regex style - RC 040721
	my ($def) = $html =~ /.*(<table class="searchtitle".*)<hr \/>/s;
	$def =~ s/([href|src]=")\//$1http\:\/\/legacy.wormbase.org\//gs;

	my $value = "<hr><h3>WormBase protein report:</h3>" . $def;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub PostHit
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $ac = $l;
	$ac =~ s/(HIT[0-9]+)\..*/$1/;
	my $cache = get_cache_root() . "/cache/$ac.hit";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,5))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			return $html;
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    	my $url= "http://www.jbirc.jbic.or.jp/hinv/spsoup/transcript_view";
	my @OutValue;
	$ac = $l;
	$ac =~ s/.+\|([A-Z][A-Z][0-9]+)\.[0-9]\|.*/$1/;
	my $ref = ['acc_id' => $ac] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(20);
	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /\<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to HIT.<BR>";
	}
	## updated to new regex style - RC 040721
	$html =~ s/.+?\<body.*?\>//si;
	$html =~ s/\<\/body.*\>.+//si;	
	$html =~ s/.+?(\<table)/$1/si;

	$html =~ s/(action\=[\"\'])(\.+)/$1http\:\/\/www.jbirc.jbic.or.jp\/hinv\/spsoup\/$2/gsi;
	$html =~ s/(src\=[\"\'])(\.+)/$1http\:\/\/www.jbirc.jbic.or.jp\/hinv\/spsoup\/$2/gsi;
	$html =~ s/(href\=[\"\'])(\.+)/$1http:\/\/www.jbirc.jbic.or.jp\/hinv\/spsoup\/$2/gsi;
	$html =~ s/tabOver\(this, \'\./tabOver\(this, 'http:\/\/www.jbirc.jbic.or.jp\/hinv\/spsoup\/./gsi;
	$html =~ s/tabOut\(this, \'\./tabOut\(this, 'http:\/\/www.jbirc.jbic.or.jp\/hinv\/spsoup\/./gsi;

	my $value = "<hr><h3>HIT protein report:</h3>" . $html;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

## returns the html from the relevant database search site or cache if it available 
## and writes the new cache file
## called by DoSearch()

sub PostScd
{
	my $html = "";
	my ($l,$protein,$refresh) = @_;
	my $cache = get_cache_root() . "/cache/$l.sgd";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,30))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			open(INPUT,"<$cache.b");
			$_ = <INPUT>;
			close(INPUT);
			if($html =~ /Moved Permanently/s)	{
				$html = "";
			}
			elsif(length($_) > 1 and not /SGD Protein Report\:/si)	{
				return $html;
			}
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			open(INPUT,"<$cache.b");
			$_ = <INPUT>;
			close(INPUT);
			if($html =~ /Moved Permanently/s)	{
				$html = "";
			}
			elsif(length($_) > 1 and not /SGD Protein Report\:/si)	{
				return $html;
			}
		}
	}
	unlink($cache . ".b");
        my $url= "http://www.yeastgenome.org/cgi-bin/locus.fpl";
	my @OutValue;
	my $ref = ['locus' => $l, 'submit' => 'Submit' , 'rm' => 'display_result'] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(20);
	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<head/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to SGD.<BR>";
	}
	## updated to new regex style - RC 040721
	$html =~ s/.+\<body.*?\>(.+)\<\/body.*/$1/si;
#	$html =~ s/.+?\<\/table.*?\>//si;
	$html =~ s/.+?(Alternative single page format)/$1/si;
	$html =~ s/.+(\<div id\=\"main\_page\"\>)/$1/si;
	$html =~ s/value=\"\/cgi-bin/value=\"http\:\/\/www.yeastgenome.org\/cgi-bin/sg;
	my $value = "<hr><h3>SGD protein report:</h3><table><tr><td>" . $html;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}


## returns the html from the relevant database search site or cache if it available 
## and writes the new cache file
## called by DoSearch()

sub PostMips
{
	my ($l,$protein,$refresh) = @_;
	my $html = "";
	my $cache = get_cache_root() . "/cache/$l.mips";
	$cache = NewCache($cache);
	my $ok = 0;
	if(not DateCheck($cache,60))	{
		$protein = 0;
		$ok = 1;
	}
	open(INPUT,"<$cache")  or $ok = 1;
	if($protein == 1){
		if(not($refresh eq "yes") and $ok == 0)	{
			local $/ = undef;  # a fancy trick - reset the record separator from \n to undef can now slurp whole file.
			$html = <INPUT>;
			close(INPUT);
			return $html;
		}
	}
	else{
		if($ok == 0)	{
			close(INPUT);
			return $html;
		}
	}
	unlink($cache . ".b");
    my $url= "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do";
	my @OutValue;
	my $ref = ['text' => $l] ;
	my $agent = LWP::UserAgent->new();
	$agent->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	my $proxy = get_proxy();
	if ($proxy=~/\w/) { $agent->proxy(['http'],$proxy); }
	$agent->agent('Mozilla/4.0 (compatible; MSIE 5.0; Windows NT; DigExt)');
	$agent->timeout(15);
	my $req = POST($url,\@OutValue,'Content_Type','form-data','Content',$ref);
	$html = $agent->request($req)->as_string();
	if(not($html =~ /<body/si))	{
		return "<BR><BR>Genomic information unavailable: cannot connect to MIPS.<BR>";
	}
	## updated to new regex style - RC 040721 - could not get it to work for all instances, so changed some back to old style
	my ($def) = $html =~ /.*(<BODY.*?)<!-- Place here your footer -->/is;
	$def =~ s/<TD.*?Search Gene\/ORF:.*?<\/TD>/<TD><\/TD>/si;
	$def =~ s/<TD.*?menu_bottom.gif.*?<\/TD>/<TD><\/TD>/si;
	$def =~ s/href=\"\//href=\"http:\/\/mips.gsf.de\//sgi;	
	$def =~ s/(href=\")([a-i|k-z])/$1http:\/\/mips.gsf.de\/genre\/proj\/yeast\/$2/sgi;	
	$def =~ s/src=\"\//src=\"http:\/\/mips.gsf.de\//sgi;
	my $value = "<hr><h3>MIPS protein report:</h3><table><tr><td>";
	$value .= $def;
	open(INPUT,">$cache");
	print INPUT $value;
	close(INPUT);
	return $value;
}

sub LoadHgnc
{
	my ($_s) = @_;
	my ($t) = $_s =~ /^(\D+)/;
	if($t eq "gi\|" or $t eq "sp\|" or $t eq "tr\|")	{
		return \%g_hgnc;
	}
	if($t eq "ENSP" and not $g_hgnc_types{"ENSP"})	{
		unless(open(IN,"<../annotation/human_hgnc.txt"))	{
			open(IN,"<human_hgnc.txt") or return \%g_hgnc;
			
		}
		$g_hgnc_types{"ENSP"} = 1;
	}
	elsif($t eq "ENSMUSP" and not $g_hgnc_types{"ENSMUSP"})	{
		unless(open(IN,"<../annotation/mouse_hgnc.txt"))	{
			open(IN,"<mouse_hgnc.txt") or return \%g_hgnc;
		}
		$g_hgnc_types{"ENSMUSP"} = 1;
	}
	elsif($t eq "AT" and not $g_hgnc_types{"ATH"})	{
		unless(open(IN,"<../annotation/ath_hgnc.txt"))	{
			open(IN,"<ath_hgnc.txt") or return \%g_hgnc;
		}
		$g_hgnc_types{"ATH"} = 1;
	}
	elsif($t eq "ENSRNOP" and not $g_hgnc_types{"ENSRNOP"})	{
		unless(open(IN,"<../annotation/rat_hgnc.txt"))	{
			open(IN,"<rat_hgnc.txt") or return \%g_hgnc;
		}
		$g_hgnc_types{"ENSRNOP"} = 1; 
	}
	elsif($t eq "FBpp" and not $g_hgnc_types{"FBpp"})	{
		unless(open(IN,"<../annotation/fly_hgnc.txt"))	{
			open(IN,"<fly_hgnc.txt") or return \%g_hgnc;
		}
		$g_hgnc_types{"FBpp"} = 1; 
	}
	elsif($t eq "ENSDARP" and not $g_hgnc_types{"ENSDARP"})	{
		unless(open(IN,"<../annotation/fish_hgnc.txt"))	{
			open(IN,"<fish_hgnc.txt") or return \%g_hgnc;
		}
		$g_hgnc_types{"ENSDARP"} = 1; 
	}
	elsif($t =~ /^Y[A-Z][LR]/ and not $g_hgnc_types{"YEAST"})	{
		unless(open(IN,"<../annotation/yeast_hgnc.txt"))	{
			open(IN,"<yeast_hgnc.txt") or return \%g_hgnc;
		}
		$g_hgnc_types{"YEAST"} = 1; 
	}
	elsif($_s =~ /[A-Z]{1,2}\d[A-Z0-9]+\.\d/ and not $g_hgnc_types{"celegans"})	{
		unless(open(IN,"<../annotation/worm_hgnc.txt"))	{
			open(IN,"<worm_hgnc.txt") or return \%g_hgnc;
		}
		$g_hgnc_types{"celegans"} = 1; 
	}
	else	{
		return \%g_hgnc;
	}
	my @values;
	my $ens;
	my @lines = <IN>;
	close(IN);
	my $l;
	foreach $l(@lines)	{
		chomp($l);
		@values = split /\t/,$l;
		$ens = @values[0];
		if(@values[1])	{
			$g_hgnc{$ens} = @values[1];
		}
		if(@values[2])	{
			$g_pubmed{$ens} = @values[2];
		}
		if(@values[3])	{
			$g_ipi{$ens} = @values[3];
		}
		if(@values[4])	{
			$g_swiss{$ens} = @values[4];
		}
	}
	return \%g_hgnc;
}

sub GetInfoBar
{
	my ($l,$gpm) = @_;
	$_ = $l;
	my $wiki = $l;
	$wiki =~ s/\|//g;
	if(/\:reversed/)	{
		return "";
	}
	my $value = $l;
	$value =~ s/\s//g;
	my $out = "";
	my $line = "";
	my $pmed = "";
	my $gpmdb = "";
	my $mrm_l = $l;
	if($mrm_l =~ /^AT[\dC]G\d/i and $mrm_l !~ /\.\d+$/)	{
		$mrm_l .= ".1";
	}
	
	$gpmdb .= " | <a GPM href=\"";
	$gpmdb .= GetHrefpSYT($mrm_l);
	$gpmdb .= "\" target=\"_info\" title=\"GPM modified peptides for $l\">pSYT</a>";
	$gpmdb .= " | <a GPM href=\"";
	$gpmdb .= GetHrefSNAP($mrm_l);
	$gpmdb .= "\" target=\"_info\" title=\"GPM amino acid polymorphism containing peptides for $l\">SNAP</a>";

	my $mrm = GetHrefMrm($mrm_l);
	my $sponsored = "";
	if($mrm)	{
		$mrm = " | <a GPM href=\"$mrm\" target=\"_info\" title=\"GPM MRM planner and proteotypic peptides for $mrm_l\">mrm</a>";
	}
	if(/ENS/ or /NEWSINFRUP/ or /GSTENP/)	{
		$out .= "| <a PROTEIN href=\"";
		$out .= GetHrefEnsembl($l);
		$out .= "\" target=\"_info\" title=\"ENSEMBL protein information for $l\">ensembl</a>";
		if(/ENSP/ or /ENSDARP/ or /ENSMUSP/ or /ENSRNOP/)	{
			$line = CheckCache($l);
			$line =~ s/\s//g;
			if(length($line) == 0 or $l =~ /ENSMUSP/ or $l =~ /ENSRNOP/)	{
				LoadHgnc($value);
				$line = $g_hgnc{$l};
				$pmed = $g_pubmed{$l};
			}
			if(length($line) > 0)	{
				if(length($out) > 0)	{
					$out .= " | ";
				}
				$out .= "<a PROTEIN href=\"";
				$out .= GetHrefNcbiType($line,"gene");
				$out .= "\" target=\"_info\" title=\"NCBI gene information for $line\">ncbi</a>";
			}
			if(length($pmed) > 0)	{
				if(length($out) > 0)	{
					$out .= " | ";
				}
				$out .= "<a GENE href=\"";
				$out .= GetHrefPubmedType($pmed,"gene");
				$out .= "\" target=\"_info\" title=\"PubMed references for the gene $line\">pubmed</a>";
			}
		}
		$_ = $l;
		if(/ENSP/)	{
			if(length($line) > 0)	{
				if(length($out) > 0)	{
					$out .= " | ";
				}
				$out .= "<a GENE href=\"";
				$out .= GetHrefNcbiType($line,"omim");
				$out .= "\" target=\"_info\" title=\"Online Mendelian Inheritance in Man for the gene $line\">omim</a> | ";
				$out .= "<a GENE href=\"";
				$out .= GetHrefNcbiType($line,"unigene");
				$out .= "\" target=\"_info\" title=\"NCBI UniGene information for the gene $line\">unigene</a> | ";
				$out .= "<a GENE href=\"";
				$out .= GetHrefNcbiType($line,"snp");
				$out .= "\" target=\"_info\" title=\"NCBI SNP information for $line\">snps</a> | ";
				$out .= "<a GENE href=\"";
				$out .= GetHrefNcbiType($line,"geo");
				$out .= "\" target=\"_info\" title=\"NCBI GEO transcriptome information for $line\">geo</a> | ";
				$out .= " <a PROTEIN href=\"";
				$out .= GetHrefProteinAtlas($_);
				$out .= "\" title=\"Human Protein Atlas tissue localization images for $_\" target=\"_info\" >hpa</a>";
				$out .= " | <a PROTEIN href=\"";
				$out .= GetHrefAntibodypedia($line);
				$out .= "\" target=\"_info\" title=\"Antibodypedia list of reagents for $line\">anti</a>";
				$out .= " | <BR>| <a PROTEIN href=\"";
				$out .= GetHrefKegg($line);
				$out .= "\" target=\"_info\" title=\"KEGG pathways containing $line\">kegg</a>";
				$out .= " | <a PROTEIN href=\"";
				$out .= GetHrefGrid($line,9606);
				$out .= "\" target=\"_info\" title=\"GRID protein-protein interactions for $line\">grid</a>";
				$out .= " | <a PROTEIN href=\"";
				$out .= GetHrefString($_);
				$out .= "\" target=\"_info\" title=\"String Interactions ($_)\">string</a>";
				$out .= " | <a PROTEIN href=\"";
				$out .= GetHrefPhosphoSite($line);
				$out .= "\" target=\"_info\" title=\"PhosphoSite protein phosphorylations for $line\">p-site</a>";
				$out .= " | <a GENE href=\"";
				$out .= GetHrefGeneCards($line);
				$out .= "\" target=\"_info\" title=\"GeneCards general information about $line\">gcard</a>";
				$out .= " | <a GENE href=\"";
				$out .= GetHrefIcgc($l);
				$out .= "\" target=\"_info\" title=\"International Cancer Genome Consortium sequence variations for $l\">icgc</a>";
				$sponsored = GetSponsored($line);
				
				my @acc = split /\t/,GetHrefUniprotType($line);
				if(@acc[1])	{
					$out .= " | <a GENE href=\"";
					$out .= "http://www.genenames.org/data/hgnc_data.php?hgnc_id=@acc[0]";
					$out .= "\" target=\"_info\" title=\"HUGO Gene Nomenclature Committee entry for $line\">hgnc</a>";
				}
			}
			$out .= " | <a PROTEIN href=\"";
			$out .= GetHrefMaxDb($_);
			$out .= "\" target=\"_info\" title=\"MaxQuant database information about $line\">maxQB</a>";
			$out .= " | <a PROTEIN href=\"";
			$out .= GetHrefPeptideAtlas($l);
			$out .= "\" target=\"_info\" title=\"PeptideAtlas database summary for $l\">pep-atlas</a>";
			$out .= " | <a PROTEIN href=\"";
			$out .= GetHrefPride($l);
			$out .= "\" target=\"_info\" title=\"PRIDE database proteomics data about $l\">pride</a>";
		}
		elsif(/ENSMUSP/)	{
			if(length($line) > 0)	{
				if(length($out) > 0)	{
					$out .= " | ";
				}
				$out .= "<a GENE href=\"";
				$out .= GetHrefNcbiType($line,"omim");
				$out .= "\" target=\"_info\" title=\"Online Mendelian Inheritance in Man entry for $line\">omim</a> | ";
				$out .= "<a PROTEIN href=\"";
				$out .= GetHrefMaxDb($_);
				$out .= "\" target=\"_info\" title=\"MaxQuant database information about $line\">maxQB</a> | ";
				$out .= "<a GENE href=\"";
				$out .= GetHrefMgi($line);
				$out .= "\" target=\"_info\" title=\"Mouse Genome Informatics entry for $line\">mgi</a> | ";
				$l = $line;
				$out .= "<a PROTEIN href=\"";
				$out .= GetHrefGrid($l,10090);
				$out .= "\" target=\"_info\" title=\"GRID protein-protein interactions for $line\">grid</a> ";
				$out .= " | <a PROTEIN href=\"";
				$out .= GetHrefAntibodypedia($line);
				$out .= "\" target=\"_info\" title=\"Antibodypedia list of reagents fro $line\">anti</a> ";
				$sponsored = GetSponsored($line);
			}
			$out .= " | <a PROTEIN href=\"";
			$out .= GetHrefPeptideAtlas($_);
			$out .= "\" target=\"_info\" title=\"PeptideAtlas database summary for $_\">pep-atlas</a>";
			$out .= " | <a PROTEIN href=\"";
			$out .= GetHrefPride($l);
			$out .= "\" target=\"_info\" title=\"PRIDE database proteomics data about $l\">pride</a>";
		}
		elsif(/ENSRNOP/)	{
		}
		$_ = $l;
	}
	if(($l =~ /^AGAP/ or $l =~ /^LOC\_Os/ or $l =~ /^BRADI/ or $l =~ /^AT[1-9C]G\d+/  or isTb($l)  or $l =~ /^Bra\d+\.\d/ 
		or /[A-Z][0-9|A-Z]+?\.[0-9]/ or /^SP[A-C][A-Z].*\-\d/ or /^PPA\d/ or /^CADAFUAP\d/) and not $l =~ /genedb/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefEnsemblMeta($l);
		$out .= "\" target=\"_info\" title=\"ENSEMBL metazoan protein entry for $l\">ensembl</a> ";
	}

	if(/^At[1-9c]g[0-9]+/i)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$l =~ s/\..+?//;
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefTair($l);
		$out .= "\" target=\"_info\" title=\"Tair protein entry for $l\">Tair</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefKegg($l);
		$out .= "\" target=\"_info\" title=\"KEGG pathways information for $l\">kegg</a> | ";
		$out .= "<a GENE href=\"";
		$out .= GetHrefNcbiType($l,"gene");
		$out .= "\" target=\"_info\" title=\"NCBI gene information for $l\">ncbi</a> | ";

		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefGrid($l,3702);
		$out .= "\" target=\"_info\" title=\"GRID protein-protein interaction information for $l\">grid</a> | ";

		$out .= "<a GENE href=\"";
		$out .= GetHrefMpss($l);
		$out .= "\" target=\"_info\" title=\"MPSS short read information for $l\">mpss</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefMipsAt($l);
		$out .= "\" target=\"_info\" title=\"Munich Information Center for Protein Sequences entry for $l\">mips</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefPride($l);
		$out .= "\" target=\"_info\" title=\"PRIDE database proteomics data about $l\">pride</a>";
	}	
	if(/ENSDARP/)	{
	}
	if(/gi\|[0-9]*\|/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefNcbi($l);
		$out .= "\" target=\"_info\" title=\"NCBI protein information for gi=$l\">ncbi</a>";
		$out .= " | <a PROTEIN href=\"";
		my ($gi) = $_ =~ /\|([0-9]+)\|/;
		$out .= "http://www.uniprot.org/mapping/?from=P_GI&to=ACC&query=$gi";
		$out .= "\" target=\"_info\" title=\"UniProt information for gi|$gi\">uniprot</a>";
		$out .= " | <a PROTEIN href=\"";
		$out .= GetHrefPride($l);
		$out .= "\" target=\"_info\" title=\"PRIDE database proteomics data about $l\">pride</a>";
	}
	if(/tgo\|.*\|/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefToxo($l);
		$out .= "\" target=\"_info\" title=\"ToxoDb protein ($l)\">toxodb</a>";
	}
	if(/sp\|[A-Z0-9_]+\|/)	{
		$out .= " | <a PROTEIN href=\"";
		$out .= GetHrefSprot($l);
		$out .= "\" target=\"_info\" title=\"SwissProt protein entry for $l\">swissprot</a>";
	}
	elsif(/tr\|[A-Z0-9_]+\|/)	{
		$out .= " | <a PROTEIN href=\"";
		$out .= GetHrefSprot($l);
		$out .= "\" target=\"_info\" title=\"SwissProt protein entry for $l\">swissprot</a>";
	}
	elsif(/plasmoDB\|/i)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefPlasmoDB($l);
		$out .= "\" target=\"_info\" title=\"PlasmoDB gene ($l)\">PlasmoDB</a>";
	}
	elsif(/DDB[0-9]+/i)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefGeneDbDdb($l);
		$out .= "\" target=\"_info\">genedb</a>";
	}	
	elsif(/genedb\|/i)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefGeneDb($_);
		$out .= "\" target=\"_info\">genedb</a>";
	}	
	elsif(/TA[0-9][0-9][0-9][0-9][0-9]/)	{
 		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefGeneDbTa($l);
		$out .= "\" target=\"_info\">genedb</a>";
	}
	elsif(/tb[0-9]+\./i)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefTigrTb($l);
		$out .= "\" target=\"_info\">tigr</a> | ";
		$out .= "<a GENE href=\"";
		$out .= GetHrefGeneDbTb($l);
		$out .= "\" target=\"_info\">genedb</a>";
	}	
	elsif(/HIT[0-9]+/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefHit($l);
		$out .= "\" target=\"_info\" title=\"HIT protein ($l)\">HIT</a>";
	}	
	elsif(/IPI[0-9]+\.?/ or /^IPI\:IPI[0-9]+\.?/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$l =~ s/^IPI\://;
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefIpi($l);
		$out .= "\" target=\"_info\" title=\"IPI protein entry for $l\">ipi</a>";
		$out .= " | <a PROTEIN href=\"";
		$out .= "http://gpmdb.thegpm.org/protein/accession/ipi\|$l";
		$out .= "\" target=\"_info\" title=\"ENSEMBL translation in GPMDB for $l\">ensembl</a>";
	}	
	elsif( /CG[0-9]+?-P[A-Z]/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefEnsembl($l);
		$out .= "\" target=\"_info\" title=\"ENSEMBL protein entry for $l\">ensembl</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefKegg($l);
		$out .= "\" target=\"_info\" title=\"KEGG pathways information for $l\">kegg</a> | ";
		$out .= "<a GENE href=\"";
		$out .= GetHrefNcbiType($l,"gene");
		$out .= "\" target=\"_info\" title=\"NCBI gene information for $l\">ncbi</a>";
	}
	if($l =~ /FBpp/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefEnsemblMeta($l);
		$out .= "\" target=\"_info\" title=\"ENSEMBL protein entry for $l\">ensembl</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefFlyBase($l);
		$out .= "\" target=\"_info\" title=\"FlyBase protein entry for $l\">flybase</a> ";
	}
	if($l =~ /DappuP/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefEnsemblMeta($l);
		$out .= "\" target=\"_info\" title=\"ENSEMBL protein entry for $l\">ensembl</a> ";
	}
	if($l =~ /[A-Z][0-9|A-Z]+?\.[0-9]/ and not (/^BRADI/ or /^AT[1-9C]G\d+/ or /genedb/))	{  # must be handled separately from the CG* labels
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefWormbase($l);
		$out .= "\" target=\"_info\" title=\"Wormbase gene entry for $l\">worm base</a>";
	}
	elsif($l =~ /^osa1/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefOsa1($l,"osa1");
		$out .= "\" target=\"_info\" title=\"OSA1 gene entry for $l\">tigr</a>";
	}
	elsif($l =~ /^hs\.\d+\|/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefHsGenome($l);
		$out .= "\" target=\"_info\" title=\"Ensembl gene location for $l\">ENSEMBL</a>";
	}
	if(/^Y[A-Z]+?[0-9]+?[A-Z]/ or /^[RQ]\d{4}/)	{
		if(length($out) > 0)	{
			$out .= " | ";
		}
		LoadHgnc($_);
		$line = $g_hgnc{$l};
		$pmed = $g_pubmed{$l};
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefScd($l);
		$out .= "\" target=\"_info\" title=\"SGD protein information for $l\">sgd</a> | ";
		if(length($pmed) > 0)	{
			$out .= "<a GENE href=\"";
			$out .= GetHrefPubmedType($pmed,"gene");
			$out .= "\" target=\"_info\" title=\"PubMed references for $line\">pubmed</a> | ";
		}
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefMips($l);
		$out .= "\" target=\"_info\" title=\"MIPS gene entry for $l\">mips</a> | ";
		$out .= "<a GENE href=\"";
		$out .= GetHrefNcbiType($l,"gene");
		$out .= "\" target=\"_info\" title=\"NCBI gene entry for $l\">ncbi</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefKegg($l);
		$out .= "\" target=\"_info\" title=\"KEGG pathways information for $l\">kegg</a> | ";
		
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefGrid($l,4932);
		$out .= "\" target=\"_info\" title=\"GRID protein-protein interactions for $l\">grid</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetHrefPride($l);
		$out .= "\" target=\"_info\" title=\"PRIDE database proteomics data about $l)\">pride</a> | ";
		
		$out .= "<a GENE href=\"";
		$out .= GetHrefNcbiType($l,"geo");
		$out .= "\" target=\"_info\" title=\"NCBI GEO transcriptome information for $l\">geo</a> | ";
		$out .= "<a PROTEIN href=\"";
		$out .= GetEorf($l);
		$out .= "\" target=\"_info\" title=\"Yeast resource center entry for $l\">yrc</a> | ";		
		$out .= "<a PROTEIN href=\"";
		$out .= "http://www.ensembl.org/Saccharomyces_cerevisiae/protview?peptide=$l&show=snps&number=on";
		$out .= "\" target=\"_info\" title=\"ENSEMBL protein entry for $l\">ensembl</a>";		
	}
	if(length($out) > 0)	{
		$out .= "$mrm | ";
	}
	$out .= "<a GPM href=\"";
	$out .= "http://$g_site{'gpmdb_url'}/thegpm-cgi/dblist_label.pl?label=$value&amp;proex=-1";
	$out .= "\" target=\"_info\" title=\"All GPMDB observations of this protein\">gpmdb</a> |";
	if($gpm =~ /GPM[0-9]+/)	{
		($gpm) = $gpm =~ /(GPM[0-9]+)/;
		$out .= " <a GPM href=\"http://$g_site{'wiki_url'}/wiki/$gpm/$wiki\" class=\"small_link\" target=\"_WIKI\" title=\"GPM wiki entry for protein $wiki in $gpm\">wiki</a>$gpmdb |"; 
	}
	else	{
		$out .= " <a GPM href=\"http://$g_site{'wiki_url'}/wiki/$wiki\" class=\"small_link\" target=\"_WIKI\" title=\"GPM wiki entry for protein accession $wiki\">wiki</a>$gpmdb |"; 
	}
	if($sponsored)	{
		$out .= "<br>Reagents: $sponsored";
	}
	return Reformat($out,$l);
}

sub Reformat
{
	my ($_in,$_label) = @_;
	$_in =~  s/\<br.*?\>//gi;
	my @links = split / \|/,$_in;
	my $l;
	my @gene;
	my @protein;
	my @other;
	my @gpm;
	foreach $l(@links)	{
		($l) = $l =~ /(<a .+?\>[^\<]+\<\/a\>)/;
		if(not $l)	{
			next;
		}
		if($l =~ /a GPM /)	{
			push(@gpm,$l);
		}
		elsif($l =~ /a PROTEIN/)	{
			push(@protein,$l);
		}
		elsif($l =~ /a GENE/)	{
			push(@gene,$l);
		}
		else	{
			push(@other,$l);
		};
	}
	my $out = qq(Click for links about <i>$_label</i>:&nbsp;&nbsp;);
	if(scalar(@protein))	{
		$out .= qq[<a class="protein" title="Show/hide links" onClick='javascript:hideBox("gene_infobar");hideBox("other_infobar");hideBox("gpm_infobar");toggleBox("protein_infobar");'>PROTEIN</a> ];
	}
	if(scalar(@gene))	{
		$out .= qq[<a class="gene" title="Show/hide links" onClick='javascript:hideBox("protein_infobar");hideBox("other_infobar");hideBox("gpm_infobar");toggleBox("gene_infobar");'>&nbsp;&nbsp;GENE&nbsp;&nbsp;</a> ];
	}
	if(scalar(@other))	{
		$out .= qq[<a class="other" title="Show/hide links" onClick='javascript:hideBox("protein_infobar");hideBox("gene_infobar");hideBox("gpm_infobar");toggleBox("other_infobar");'>OTHER</a> ];
	}
	if(scalar(@gpm))	{
		$out .= qq[<a class="gpm" title="Show/hide links" onClick='javascript:hideBox("protein_infobar");hideBox("gene_infobar");hideBox("other_infobar");toggleBox("gpm_infobar");'>&nbsp;GPMDB&nbsp;</a> ];
	}
	if(scalar(@protein))	{
		$out .= qq(<div id="protein_infobar" class="links" style="display: none; background-color: #D6F8DE; width: 605px; border: 1px solid #669966;">);
		my $row = 0;
		$out .= qq(<table width="600" cellspacing="2" cellpadding="1" border="0"><tr><td colspan="3">Links to protein-specific external resources</td>
		<td align="right" valign="middle"><img src="/llama-magic-html/ffa_expanded.gif" onClick="javascript: hideBox('protein_infobar');" border="0" title="hide links" /></td></tr><tr>);
		my $text;
		foreach $l(@protein)	{
			($text) = $l =~ /title\=\"(.+?)\"/;
			if($row == 0)	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			elsif(not $row % 2)	{
				$out .= qq(</tr><tr><td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			else	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			$row++;
		}
		if(not $row % 2)	{
			$out .= qq(</tr>);
		}
		else	{
			$out .= qq(<td width="75"></td><td width="213"></td></tr>);
		}
		$out .= qq(</table></div>);
	}
	if(scalar(@gene))	{
		$out .= qq(<div id="gene_infobar" class="links" style="display: none; background-color: #F8E9FC; width: 605px; border: 1px solid #669966;">);
		my $row = 0;
		$out .= qq(<table width="600" cellspacing="2" cellpadding="1" border="0"><tr><td colspan="3">Links to gene-specific external resources</td>
		<td align="right" valign="middle"><img src="/llama-magic-html/ffa_expanded.gif" onClick="javascript: hideBox('gene_infobar');" border="0" title="hide links" /></td></tr><tr>);
		my $text;
		foreach $l(@gene)	{
			($text) = $l =~ /title\=\"(.+?)\"/;
			if($row == 0)	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			elsif(not $row % 2)	{
				$out .= qq(</tr><tr><td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			else	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			$row++;
		}
		if(not $row % 2)	{
			$out .= qq(</tr>);
		}
		else	{
			$out .= qq(<td width="75"></td><td width="213"></td></tr>);
		}
		$out .= qq(</table></div>);
	}
	if(scalar(@other))	{
		$out .= qq(<div id="other_infobar" class="links" style="display: none; background-color: #FFEAEA; width: 605px; border: 1px solid #669966;">);
		my $row = 0;
		$out .= qq(<table width="600" cellspacing="2" cellpadding="1" border="0"><tr><td colspan="3">Links to other external resources</td>
		<td align="right" valign="middle"><img src="/llama-magic-html/ffa_expanded.gif" onClick="javascript: hideBox('other_infobar');" border="0" title="hide links" /></td></tr><tr>);
		my $text;
		foreach $l(@other)	{
			($text) = $l =~ /title\=\"(.+?)\"/;
			if($row == 0)	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			elsif(not $row % 2)	{
				$out .= qq(</tr><tr><td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			else	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			$row++;
		}
		if(not $row % 2)	{
			$out .= qq(</tr>);
		}
		else	{
			$out .= qq(<td width="75"></td><td width="213"></td></tr>);
		}
		$out .= qq(</table></div>);
	}
	if(scalar(@gpm))	{
		$out .= qq(<div id="gpm_infobar" class="links" style="display: none; background-color: #FFF7B7; width: 605px; border: 1px solid #669966;">);
		my $row = 0;
		$out .= qq(<table width="600" cellspacing="2" cellpadding="1" border="0"><tr><td colspan="3">Links to GPMDB internal resource</td>
		<td align="right" valign="middle"><img src="/llama-magic-html/ffa_expanded.gif" onClick="javascript: hideBox('gpm_infobar');" border="0" title="hide links" /></td></tr><tr>);
		my $text;
		foreach $l(@gpm)	{
			($text) = $l =~ /title\=\"(.+?)\"/;
			if($row == 0)	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			elsif(not $row % 2)	{
				$out .= qq(</tr><tr><td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			else	{
				$out .= qq(<td valign="top" align="right" width="75">$l:</td><td valign="top" align="left" width="213">$text</td>);
			}
			$row++;
		}
		if(not $row % 2)	{
			$out .= qq(</tr>);
		}
		else	{
			$out .= qq(<td width="75"></td><td width="213"></td></tr>);
		}
		$out .= qq(</table></div>);
	}
	return $out;		
}
	
sub GetHrefIcgc
{
	my ($_acc) = @_;
	return "http://dcc.icgc.org/martsearch/#!/?q=$_acc";
}

sub GetHrefMaxDb
{
	my ($l) = @_;
	return "http://www.biochem.mpg.de/maxqb/mxdb/protein/show?sourceId=$l";
}
			
sub GetHrefMrm
{
	my ($l,$t) = @_;
	if($l =~ /gi\|/)	{
		return "";
	}
	my ($p) = $l =~ /(^[A-Z]+)/i;
	my $ua = LWP::UserAgent->new;
	$ua->cookie_jar(HTTP::Cookies->new(file => $g_cookie_jar,autosave => 1));

	$ua->timeout(5);
	my $uau = $g_site{'gpmdb_user'};
	my $uap = $g_site{'gpmdb_pwd'};
	my $urla = "http://";
	my $header;
	my $check_url = "http://mrm.thegpm.org/thegpm-cgi/ajax_server_mrm.pl?target=label_check&data=$l";
	if($uau and $uap)	{
		$urla .= "$uau:$uap@";
	}
	if($t)	{
		return "http://$g_site{'mrm_url'}/thegpm-cgi/peak_search.pl?pmass=&pmrange=0.1&pattern=$p&label=$l&seq=&unique=on&fmrange=0.1&sort=int&submit=Search";
	}		
	$urla .= get_server_name() . "/thegpm-cgi/request_server_mrm.pl?target=label_check&data=$l";
	$header = $ua->head($check_url);

	if ($header->{_rc} != 403) {  # label found in MRM data

		return "http://$g_site{'mrm_url'}/thegpm-cgi/peak_search.pl?pmass=&pmrange=0.1&pattern=$p&label=$l&seq=&unique=on&fmrange=0.1&sort=int&submit=Search";
	}
	return "";
}

sub GetHrefUniprot
{
	my ($l) = @_;
	return "http://www.uniprot.org/uniprot/?query=$l";
}

sub GetHrefNextprot
{
	my ($l) = @_;
	return "http://www.nextprot.org/db/search#$l";
}

sub GetHrefpSYT
{
	my ($l) = @_;
	return "http://$g_site{'psyt_url'}/thegpm-cgi/dblist_pep_modmass.pl?label=$l&modmass=80\@STY&display=0";
}

sub GetHrefSNAP
{
	my ($l) = @_;
	return "http://$g_site{'snap_url'}/thegpm-cgi/dblist_protein_mut.pl?label=$l";
}

sub GetHrefPride
{
	my ($l) = @_;
	$l =~ s/\|$//;
	return "http://www.ebi.ac.uk/pride/searchSummary.do?queryTypeSelected=identification accession number&identificationAccessionNumber=$l";
}

sub GetHrefProteinAtlas
{
	my ($l) = @_;
	return "http://www.proteinatlas.org/search/$l";
}

sub GetHrefPeptideAtlas
{
	my ($l) = @_;
	if($l =~ /^ENSP\d/)	{
		return "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=335&protein_name=$l&action=QUERY";
	}
	if($l =~ /^ENSMUSP\d/)	{
		return "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=317&protein_name=$l&action=QUERY";
	}
	if($l =~ /^Y[A-Z][LR]/)	{
		return "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=180&protein_name=$l&action=QUERY";
	}
	return "";
}

sub GetHrefPlasmoDB
{
	my ($l) = @_;
	$l =~ s/plasmodb\|//i;
	return "http://www.plasmodb.org/plasmo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&primary_key=$l";
}

sub GetHrefToxo
{
	my ($l) = @_;
	$l =~ s/tgo\|(.+?)\|/$1/i;
	return "http://www.toxodb.org/toxo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&primary_key=$l&project_id=ToxoDB";
}


sub GetHrefHapMap
{
	my ($l) = @_;
	return "http://www.hapmap.org/cgi-perl/gbrowse/hapmap/?pluggin=PairplotAnnotator;name=$l";
}

sub GetHrefEnsembl
{
	my ($l) = @_;
    ##added 20050215
    my ($label) = $l =~ /(ENS\w*[PG])\d+/; 
    if($label){
		return "http://www.ensembl.org/$g_ens{$label}/protview?peptide=$l&show=snps&number=on";
	}
	elsif($l =~ /CG[0-9]+?-P[A-Z]/)	{
		return "http://www.ensembl.org/Drosophila_melanogaster/protview?peptide=$l&show=snps&number=on";
	}
	elsif($l =~ /FBpp/)	{
		return "http://www.ensembl.org/Drosophila_melanogaster/protview?peptide=$l&show=snps&number=on";
	}
	elsif($l =~ /[A-Z][0-9|A-Z]+?\.[0-9]/)	{
		return "http://www.ensembl.org/Caenorhabditis_elegans/protview?peptide=$l&show=snps&number=on";
	}
	return "http://www.ensembl.org/Homo_sapiens/protview?peptide=$l&show=snps&number=on";
}

sub GetHrefEnsemblMeta
{
	my ($l) = @_;
    ##added 20050215
	my ($label) = $l =~ /([A-Z\_]+)\d/i;
	if($label =~ /^LOC\_Os/  or $label =~ /^AT/ or $label =~ /BRADI/){
		$l =~ s/\-P$/-TAIR/;
		return "http://plants.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;t=$l";
	}
	elsif($label =~ /^AA[QZ]/ or $label =~ /^EAN/ or $label =~ /^CAJ/){
		return "http://protists.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;t=$l";
	}
	elsif($label =~ /^Bra/){
		return "http://plants.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;t=$l";
	}
	elsif($l =~ /^SP[A-C][A-Z].*\-\d/)	{
		return "http://fungi.ensembl.org/Schizosaccharomyces_pombe/Transcript/ProteinSummary?db=core;t=$l";
	}	
	elsif($l =~ /^CADAFUAP\d/)	{
		return "http://fungi.ensembl.org/Aspergillus_fumigatus/Transcript/ProteinSummary?db=core;t=$l";
	}	
	elsif($l =~ /[A-Z][0-9|A-Z]+?\.[0-9]/)	{
		return "http://metazoa.ensembl.org/Caenorhabditis_elegans/Transcript/ProteinSummary?db=core;t=$l";
	}
	return "http://metazoa.ensembl.org/$g_meta{$label}/Transcript/ProteinSummary?db=core;t=$l";
}

sub GetHrefScd
{
	my ($l) = @_;
	return "http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=$l"; 
}

sub GetHrefKegg
{
	my ($l) = @_;
	return "http://www.genome.ad.jp/dbget-bin/www_bfind_sub?max_hit=1000&dbkey=genes&mode=bfind&keywords=$l";
}

sub GetHrefAntibodypedia	{
	my ($l) = @_;
	return "http://www.antibodypedia.com/genes.php?query=$l";
}

sub GetHrefHmdb
{
	my ($l) = @_;
	return "http://www.hmdb.ca/search/search?query=$l";
}

sub GetHrefFlyBase
{
	my ($l) = @_;
	return "http://flybase.org/reports/$l.html";
}
sub GetHrefHprd
{
	my ($l) = @_;
	return "http://www.hprd.org/resultsQuery?multiplefound=&prot_name=$l&external=Ref_seq&accession_id=&hprd=&gene_symbol=&chromo_locus=&function=&ptm_type=&localization=&domain=&motif=&expression=&prot_start=&prot_end=&limit=0&mole_start=&mole_end=&disease=&query_submit=Search";
}

sub GetEorf
{
	my ($l) = @_;
	return "http://www.yeastrc.org/pdr/yeastProteinRedirect.do?acc=$l";
}

sub GetHrefNcbi
{
	my ($l) = @_;
	my $gi = $l;
	my $nuc = 0;
	if($gi =~ /gb\|/)	{
		$nuc = 1;
	}
	$gi =~ s/gi\|//s;
	$gi =~ s/\|.*//s;
	if($nuc == 1)	{
		return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=Nucleotide&val=$gi";
	}
	return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=$gi";
}

sub GetHrefSprot
{
	my ($l) = @_;
	my $gi = $l;
	if($gi =~ /^sp\|/)	{
		$gi =~ s/sp\|//s;
	}
	elsif($gi =~ /^tr\|/)	{
		$gi =~ s/tr\|//s;
	}
	$gi =~ s/\|.*//s;
	return "http://www.uniprot.org/uniprot/$gi.txt";
}

sub GetHrefMgi
{
	my ($l) = @_;
	return "http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=searchTool&query=$l&selectedQuery=Genes+and+Markers";
}

sub GetHrefRgd
{
	my ($l) = @_;
	return "http://rgd.mcw.edu/generalSearch/RgdSearch.jsp?searchKeyword=$l&quickSearch=1";
}

sub GetHrefNcbiType
{
	my ($l,$t) = @_;
	$l =~ s/:/\%20/g;
	if($t eq "omim")	{
		return "http://www.omim.org/search?search=$l";
	}
	return "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=$t&term=$l";
}

sub GetHrefPubmedType
{
	my ($l,$t) = @_;
	return "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Link&LinkName=gene_pubmed&from_uid=$l";
}

sub GetHrefUniprotType
{
	my ($l) = @_;
	my @acc;
	my $found = 0;
	if(-e "gdl.txt")	{
		open(FILE,"<gdl.txt") or return "";
	}
	else	{
		open(FILE,"<../annotation/gdl.txt") or return "";
	}
	while(<FILE>)	{
		if(/\t$l\t/)	{
			$found = 1;
			last;
		}
	}
	close(FILE);
	if($found)	{
		return $_;
	}
	return "";
}

sub GetHrefOsa1
{
	my ($l,$t) = @_;
	$l =~ s/osa1//;
	return "http://www.tigr.org/tigr-scripts/euk_manatee/shared/ORF_infopage.cgi?db=$t&orf=$l";
}

sub GetHrefHsGenome
{
	my ($l) = @_;
	my ($chr) = $l =~ /hs\.(\d+)\|/;
	my ($start) = $l =~ /hs\.\d+\|(\d+)\-\d+\|/;
	my ($end) = $l =~ /hs\.\d+\|\d+\-(\d+)\|/;
	return "http://www.ensembl.org/Homo_sapiens/contigview?chr=$chr&region=&start=$start&end=$end";
}

sub GetHrefHit
{
	my ($l,$t) = @_;
	$l =~ s/.+\|([A-Z][A-Z][0-9]+)\.[0-9]\|.*/$1/;
	return "http://www.jbirc.jbic.or.jp/hinv/soup/pub_Detail.pl?acc_id=$l";
}

sub GetHrefIpi
{
	my ($l,$t) = @_;
	if($l =~ /^IPI\:/)	{
		$l =~ s/\|.*//;
		return "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-e+[$l]+-vn+2";
	}
	$l =~ s/\..*//;
	return "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-e+[IPI-acc:$l]+-vn+2";
}

sub GetHrefTair
{
	my ($l) = @_;
	return "http://arabidopsis.org/servlets/TairObject?type=locus&name=$l";
}

sub GetHrefMipsAt
{
	my ($l) = @_;
	return "http://mips.helmholtz-muenchen.de/plant/athal/searchjsp/index.jsp?searchge:text1=$l&searchge:searchby=ByFreetext&searchge:cb1=submit&searchge=searchge";
}
sub GetHrefBind
{
	my ($l) = @_;
	return "http://bind.ca/Action?textquery=RecordType: (interaction complex pathway )+AND+$l";
}
sub GetHrefPhosphoSite{

	my ($l) = @_;
	return "http://www.phosphosite.org/simpleSearchSubmitAction.do?queryId=-1&from=0&searchStr=$l&x=0&y=0";
}
sub GetHrefGeneCards{

	my ($l) = @_;
	return "http://www.genecards.org/cgi-bin/carddisp.pl?gene=$l";
}

sub GetHrefGrid
{
	my ($l,$i) = @_;
	return "http://www.thebiogrid.org/search.php?keywords=$l&organismid=$i";
}

sub GetHrefString
{
	my ($l) = @_;
	return "http://string-db.org/newstring_cgi/show_network_section.pl?identifier=$l";
}

sub GetHrefMpss
{
	my ($l) = @_;
	return "http://mpss.udel.edu/at/GeneAnalysis.php?featureName=$l";
}

sub GetHrefTigrAt
{
	my ($l) = @_;
	return "http://www.tigr.org/tigr-scripts/euk_manatee/shared/ORF_infopage.cgi?db=ath1&orf=$l";
}

sub GetHrefTigrTb
{
	my ($l) = @_;
	return "http://www.tigr.org/tigr-scripts/euk_manatee/shared/ORF_infopage.cgi?db=tba1&orf=$l";
}

sub GetHrefGeneDb
{
	my ($l) = @_;
	($l) = $l =~ /\|(.+?)\|/;
	return "http://www.genedb.org/gene/$l";
}

sub GetHrefGeneDbTb
{
	my ($l) = @_;
	return "http://www.genedb.org/genedb/Search?submit=Search+for&name=$l&organism=tryp&desc=yes&wildcard=yes";
}

sub GetHrefGeneDbDdb
{
	my ($l) = @_;
	return "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=dicty&name=$l&submit=Search";
}

sub GetHrefGeneDbSp
{
	my ($l) = @_;
	return "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=pombe&name=$l&submit=Search";
}

sub GetHrefGeneDbTa
{
	my ($l) = @_;
   	return "http://www.genedb.org/genedb/Search?submit=Search+for&name=$l&organism=annulata&desc=yes&wildcard=yes";
}

sub GetHrefMips
{
	my ($l) = @_;
	return "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=$l";
}

sub GetHrefWormbase
{
	my ($l) = @_;
	return "http://legacy.wormbase.org/db/gene/gene?name=$l";
}

sub CheckCache
{
	my ($l) = @_;	
	my $cache = get_cache_root() . "/cache/$l.ensembl";
	$cache = NewCache($cache);
	my $html = "";
	my $ok = 0;
	open(INPUT,"<$cache")  or return $html;
	while(<INPUT>)	{
		if(/\<th.*?\>Peptide\<\/th\>/)	{
			$_ = <INPUT>;
			if(/HUGO ID/i or /ZFIN\_ID ID/i)	{
				s/.*?\<b\>(.+)\<\/b\>.*/$1/i;
				$html = $_;
				close(INPUT);
				return $html;
			}
		}
		if(/HUGO ID/i or /ZFIN\_ID ID/i)	{
			s/.*?\<strong>(.+)\<\/strong\>.*/$1/i;
			$html = $_;
			close(INPUT);
			return $html;
		}
		if(/HGNC Symbol ID/i)	{
			($html) = /\<strong\>(.+)\<\/strong\>/i; 
			close(INPUT);
			return $html;
		}
		if(/MGI:[0-9]+/)	{
			($html) = /(MGI:[0-9]+)/; 
			close(INPUT);
			return $html;
		}
	}
	close(INPUT);
	return $html;
}	

sub GetMarkup
{
	my ($l) = @_;	
	if($l =~ /\:reversed/)	{
		return "";
	}
	my $cache = get_cache_root() . "/cache/$l.ensembl";
	$cache = NewCache($cache);
	my $html = "";
	my $ok = 0;
	open(INPUT,"<$cache")  or return $html;
	while(<INPUT>)	{
		if(/\<pre.*?\>/ and not /error/)	{
			if(not /span style\=\"color\:/)	{
				$_ = <INPUT>;
			}
			else	{
				s/.*\<pre.+?\>//;
			}
			while(not/\<\/pre/)	{
				if(/span/)	{
					s/background\-color\:ff9999/font-style:normal;text-decoration:underline;text-transform:uppercase;cursor:pointer;background-color:#ffaaaa/g;
					s/background\-color\:99ff99/font-style:normal;text-decoration:overline;text-transform:uppercase;cursor:pointer;background-color:#aaffaa/g;
					s/          / /;
					$html .= $_;
				}
				$_ = <INPUT>;
			}
			close(INPUT);
			return $html;
		}
	}
	close(INPUT);
	return $html;
}	

sub cover_svg
{
	my ($e,$p,$start,$end,$sc,$exp,$u,$base_color,$_un,$_height) = @_;
	if(length($base_color) == 0)	{
		$base_color = "red";
	}
	my $height = 15;
	if($_height)	{
		$height = $_height;
	}
	my $s = 1;

open(OUTPUT,">$p");
print OUTPUT <<End_of_svg;
<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 20010904//EN"
"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="$sc" height="$height">
	<defs>
		<filter id="MyFilter" filterUnits="userSpaceOnUse">
			<!--Copyright 1999 Adobe Systems. You may copy, modify, and distribute this file, if you include this notice & do not charge for the distribution. This file is provided "AS-IS" without warranties of any kind, including any implied warranties.-->
			<feGaussianBlur in="SourceAlpha" stdDeviation="3" result="blur"/>
			<feOffset in="blur" dx="3" dy="2" result="offsetBlurredAlpha"/>
			<feSpecularLighting in="blur" surfaceScale="3" specularConstant="0.6" specularExponent="20" lightColor="rgb(128,128,128)" result="specularOut">
				<feDistantLight azimuth="235" elevation="40"/>
			</feSpecularLighting>
			<feComposite in="specularOut" in2="SourceAlpha" operator="in" result="specularOut"/>
			<feComposite in="SourceGraphic" in2="specularOut" operator="arithmetic" k1="0" k2="1" k3="1" k4="0" result="litPaint"/>
			<feMerge>
				<feMergeNode in="offsetBlurredAlpha"/>
				<feMergeNode in="litPaint"/>
			</feMerge>
		</filter>
	</defs>

End_of_svg

	my $a = 0;
	my $xoff = 0;
	my $y = int($height/2 +0.5);
	if($sc == 0)	{
		$sc = 600;
	}
	my $xlength = $sc-$xoff;
	my $scale = $xlength/$e;
	my $x1 = int($xoff + $scale*($s-1) + 0.5);
	my $x2 = int($xoff + $scale*($e-1) + 0.5);
	my %pairs;
#
# make a single enclosing <a> tag for the whole graphic, including a terrible kludge to make IE behave.
#

	my $browser = $ENV{'HTTP_USER_AGENT'};
#print  qq(<script>alert('browser is $browser');</script>);

	if ($browser =~ /msie/ig and length($u) > 0) { # explorer, so write tag for it

		print OUTPUT "<a xlink:href=\"$u\" xlink:show=\"new\" xlink:target=\"_blank\">\n";

	} elsif(length($u) > 0) { # write tag for everything else

		print OUTPUT "<a xlink:href=\"$u\" xlink:show=\"new\" target=\"_blank\" xlink:target=\"_blank\">\n";

	} # end if
	my $value = "<line x1=\"$x1\" y1=\"$y\" x2=\"$x2\" y2=\"$y\" style=\"stroke:black; stroke-width:1;\"/>\n";

	my $cent = 100;
	my $ycent1 = $y+3;
	my $ycent2 = $y-3;
	while($cent < $e)	{
		$x2 =  int($xoff + $scale*($s-1+$cent) + 0.5);
		if($cent % 500)	{
			$ycent1 = $y+3;
			$ycent2 = $y-3;
			print OUTPUT "<line x1=\"$x2\" y1=\"$ycent1\" x2=\"$x2\" y2=\"$ycent2\" style=\"stroke:black; stroke-width:1;\"/>";
		}
		else	{
			$ycent1 = $y+3;
			$ycent2 = $y-3;
			print OUTPUT "<line x1=\"$x2\" y1=\"$ycent1\" x2=\"$x2\" y2=\"$ycent2\" style=\"stroke:black; stroke-width:2;\"/>";
		}
		$cent += 100;
	}
	print OUTPUT $value;
	my $op;
	my $color;
	my $convert = -1.0/(8.0*log(10));
	my $cx;
	my $cy;
	my $rx;
	my $ry;
	my $width = int($height/2 + 1);
	print OUTPUT "<g style=\"stroke-width:$width; stroke:green;\">\n";
	while($a < $e)	{
		($value,$op) = split / /,$$_un{$a};
		my $st = $a+1;
		my $en;
		if($value)	{
			$x1 = $xoff + $scale*($a);
			$x1 = int($x1 + 0.5);
			if($a + $value < $e)	{
				$x2 = $xoff + $scale*($a+$value);
				$en = $a+$value;
			}
			else	{
				$x2 = $xoff + $scale*($e-1);
				$en = $e;
			}
			$x2 = int($x2 + 0.5);
			if($x2 - $x1 < 1)	{
				$x2 = $x1 + 1;
			}
			if(length($u) > 0)	{
#				print OUTPUT "<a xlink:href=\"$u\" xlink:show=\"new\">";
			}
			if(not $op)	{	
				print OUTPUT "<line xlink:title=\"$st-$en\" style=\"opacity:0.3; stroke-width:$width; stroke:green;\" x1=\"$x1\" x2=\"$x2\" y1=\"$y\" y2=\"$y\"><title>$st-$en</title></line>";
			}
			else	{
				if($op < 0.29)	{
					print OUTPUT "<line xlink:title=\"$st-$en\" style=\"opacity:0.3; stroke-width:$width; stroke:cyan;\" x1=\"$x1\" x2=\"$x2\" y1=\"$y\" y2=\"$y\"><title>$st-$en</title></line>";
				}
				else	{
					print OUTPUT "<line xlink:title=\"$st-$en\" style=\"opacity:$op; stroke-width:$width; stroke:green;\" x1=\"$x1\" x2=\"$x2\" y1=\"$y\" y2=\"$y\"><title>$st-$en</title></line>";
				}
			}
			if(length($u) > 0)	{
#				print OUTPUT "</a>";
			}
			print OUTPUT "\n";
			$a += $value;
		}
		else	{
			$a++;
		}
	}
	print OUTPUT "</g>\n";

	$a = 0;
	print OUTPUT "<g style=\"stroke-width:$width; stroke:$base_color;\">\n";
	my @pro_len = ((0) x $e);
	my $res_ctr;

	foreach $value(@$start)	{

		if($pairs{$value} != @$end[$a])	{

			for ($res_ctr=$value; $res_ctr<=@$end[$a]; $res_ctr++) {  # set covered section

				$pro_len[$res_ctr] = 1;

			}  # end for

			$x1 = $xoff + $scale*($value-1);
			$x1 = int($x1 + 0.5);
			$x2 = $xoff + $scale*(@$end[$a]-1);
			$x2 = int($x2 + 0.5);
			if($x2 - $x1 < 1)	{
				$x2 = $x1 + 1;
			}
			$pairs{$value} = @$end[$a];
			my $st = $value;
			my $en = @$end[$a];	
			if(not $exp)	{	
				if(length($u) > 0)	{
#					print OUTPUT "<a xlink:href=\"$u\" xlink:show=\"new\">";
				}
				print OUTPUT "<line x1=\"$x1\" x2=\"$x2\" y1=\"$y\" y2=\"$y\"/>";
				if(length($u) > 0)	{
#					print OUTPUT "</a>";
				}
				print OUTPUT "\n";
			}
			else	{
				if(@$exp[$a] > 0){
					$color = $convert*log(@$exp[$a]);
				}
				else{
					$color = 1;
				}
				if($color > 1)	{
					$color = 1;
				}
				if($color < .1)	{
					$color = .1;
				}
				$color = sprintf("%.2f",$color);
				if(length($u) > 0)	{
#					print OUTPUT "<a xlink:href=\"$u\" xlink:show=\"new\">";
				}
				print OUTPUT "<line xlink:title=\"$st-$en\" style=\"opacity:$color; stroke-width:$width; stroke:$base_color;\" x1=\"$x1\" x2=\"$x2\" y1=\"$y\" y2=\"$y\"><title>$st-$en</title></line>\n";
				if(length($u) > 0)	{
#					print OUTPUT "</a>";
				}
				print OUTPUT "\n";
			}	
		}
		$a++;
	}

	$res_ctr = 0;

	foreach (@pro_len) {  # count up covered residues

		if ($_ == 1) { 

			$res_ctr++;

		}  # end if

	}  # end foreach

	my $pro_coverage = sprintf("%.1f", ($res_ctr/$e)*100);

#	print OUTPUT "</g>\n</svg>\n";
	if(length($u) > 0)	{
		print OUTPUT "</g>\n</a>\n</svg>\n";
	}
	else	{
		print OUTPUT "</g>\n</svg>\n";
	}
	close OUTPUT;

	return $pro_coverage;

}

sub print_pgo
{
	my ($pro,$url,$proex,$npep,$lt) = @_;
	my $a = 0;
	my $r;
	foreach $a(@$pro)	{
		if($a =~ /^ENSP/)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=human&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a>&nbsp;|&nbsp;";
			$r .= "<a href=\"/thegpm-cgi/ptissue.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=human&amp;ltype=$lt\" class=\"small_link\" title=\"BRENDA classification diagram\">BTO</a>&nbsp;|&nbsp;";
			$r .= "<a href=\"/thegpm-cgi/ppath.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=human&amp;ltype=$lt\" class=\"small_link\" title=\"KEGG pathways analysis\">path</a>&nbsp;|&nbsp;";
			$r .= "<a href=\"/thegpm-cgi/pppi.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=human&amp;ltype=$lt\" class=\"small_link\" title=\"Protein-protein interaction analysis\">ppi</a>&nbsp;|&nbsp;";
			$r .= "<a href=\"/thegpm-cgi/pdomain.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=human&amp;ltype=$lt\" class=\"small_link\" title=\"Protein domain analysis\">doms</a>&nbsp;|&nbsp;";
			$r .= "<a href=\"/thegpm-cgi/psap.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=human&amp;ltype=$lt\" class=\"small_link\" title=\"Single amino acid polymorphisms\">snaps</a> | ";
			last;
		}
		if($a =~ /^FBpp/)	{
			$r = "<a href=\"/thegpm-cgi/psap.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=human&amp;ltype=$lt\" class=\"small_link\" title=\"Single amino acid polymorphisms\">snaps</a> | ";
			last;
		}
		if($a =~ /^Y[A-Z][A-Z]/ or $a =~ /^[RQ]\d{4}/)	{
			$r =  "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=yeast&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | <a href=\"/thegpm-cgi/ppath.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=yeast\" class=\"small_link\" title=\"KEGG pathways analysis\">path</a> | ";
			$r .= "<a href=\"/thegpm-cgi/pppi.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=yeast&amp;ltype=$lt\" class=\"small_link\" title=\"Protein-protein interaction analysis\">ppi</a>&nbsp;|&nbsp;";
			$r .= "<a href=\"/thegpm-cgi/pdomain.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=yeast&amp;ltype=$lt\" class=\"small_link\" title=\"Protein domain analysis\">doms</a>&nbsp;|&nbsp;";
			$r .= "<a href=\"/thegpm-cgi/psap.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=yeast&amp;ltype=$lt\" class=\"small_link\" title=\"Single amino acid polymorphisms\">snaps</a> | ";
			last;
		}
		if($a =~ /^ENSMUSP/)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=mouse&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | ";
			$r .= "<a href=\"/thegpm-cgi/ppath.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=mouse&amp;ltype=$lt\" class=\"small_link\" title=\"KEGG pathways analysis\">path</a> | ";
			$r .= "<a href=\"/thegpm-cgi/pdomain.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=mouse&amp;ltype=$lt\" class=\"small_link\" title=\"Protein domain analysis\">doms</a> | ";
			$r .= "<a href=\"/thegpm-cgi/psap.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=mouse&amp;ltype=$lt\" class=\"small_link\" title=\"Single amino acid polymorphisms\">snaps</a> | ";
			last;
		}
		if($a =~ /^ENSBTAP/)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=cow&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | ";
			last;
		}
		if($a =~ /^ENSDARP/)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=fish&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | ";
			$r .= "<a href=\"/thegpm-cgi/ppath.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=fish&amp;ltype=$lt\" class=\"small_link\" title=\"KEGG pathways analysis\">path</a> | ";
			last;
		}
		if($a =~ /^ENSRNOP/)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=rat&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | ";
			$r .= "<a href=\"/thegpm-cgi/ppath.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=rat&amp;ltype=$lt\" class=\"small_link\" title=\"KEGG pathways analysis\">path</a> | ";
			$r .= "<a href=\"/thegpm-cgi/pdomain.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=rat&amp;ltype=$lt\" class=\"small_link\" title=\"Protein domain analysis\">doms</a> | ";
			$r .= "<a href=\"/thegpm-cgi/psap.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=rat&amp;ltype=$lt\" class=\"small_link\" title=\"Single amino acid polymorphisms\">snaps</a> | ";
			last;
		}
		if($a =~ /^ENSCAFP/)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=dog&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | ";
			last;
		}
		if($a =~ /^AT.G/i)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=ath&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | ";
			$r .= "<a href=\"/thegpm-cgi/ppath.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=ath&amp;ltype=$lt\" class=\"small_link\" title=\"KEGG pathways analysis\">path</a> | ";
			$r .= "<a href=\"/thegpm-cgi/psap.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=ath&amp;ltype=$lt\" class=\"small_link\" title=\"Single amino acid polymorphisms\">snaps</a> | ";
			last;
		}
		if($a =~ /^[A-Z][0-9|A-Z]+?\.[0-9]/)	{
			$r = "<a href=\"/thegpm-cgi/pgo.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=worm&amp;ltype=$lt\" class=\"small_link\" title=\"GO classification diagram\">GO</a> | ";
			$r .= "<a href=\"/thegpm-cgi/ppath.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;taxa=worm&amp;ltype=$lt\" class=\"small_link\" title=\"KEGG pathways analysis\">path</a> | ";
			last;
		}
	}
	$r .= "<a href=\"/thegpm-cgi/pepdelta.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;style=1&amp;ltype=$lt\" class=\"small_link\" title=\"Parent ion mass error distribution\">mh</a> | ";
	$r .= "<a href=\"/thegpm-cgi/pepdelta.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;style=3&amp;ltype=$lt\" class=\"small_link\" title=\"Ion zeta distribution\">&zeta;</a> | ";
	return $r;
}

sub biopixie_form
{
	my ($_p) = @_;
	my $proteins;
	my $a;
	foreach $a(@$_p)	{
		if($a =~ /^Y[A-Z][LR]\d/ or $a =~ /^[RQ]\d{4}/)	{
			$proteins .= "$a,";
		}
		if(length($proteins) > 990)	{
			last;
		}
	}
	print qq(<div style="display:none"> <form name="biopixie" action="HTTP://pixie.princeton.edu/pixie/graph.php" enctype="multipart/form-data" method="post">\n);
	print qq(<input type="hidden" name="graph_genes" value="$proteins">\n);
	print qq(</form></div>\n);
}

sub simple_model_toolbar {
#
# This function takes a single GPM number as an argument and returns a string consisting of HTML to describe
#	the links for this model that would be found on the top of the plist.pl page.
#
# Inputs: a GPM number.
#
# Outputs: text describing HTML
#

	my $gpmnum = shift; # read in from argument
	my $result; # scalar to hold text
	my $full_gpm_path = '/gpm/archive/' . substr($gpmnum, 3, 3) . "/$gpmnum" . '.xml';

	$result .= qq(<a href="/thegpm-cgi/plist.pl?npep=0&path=$full_gpm_path&amp;proex=-1" class="small_link" title="Summary navigation page for this model">model</a> | 
<a href="/thegpm-cgi/pintersect.pl?path=$full_gpm_path&amp;proex=-1&amp;npep=0" title="Locate all models containing a set of proteins">context</a> | 
<a href="/thegpm-cgi/pgroup.pl?path=$full_gpm_path&amp;proex=-1&amp;npep=0" title="Show protein grouping based on peptide homology">group</a> | 
<a href="/thegpm-cgi/pgel.pl?path=$full_gpm_path&amp;proex=-1&amp;npep=0" title="A 1D and 2D PAGE simulation of this model">gel</a> | 
<a href="/thegpm-cgi/pchip.pl?path=$full_gpm_path&amp;proex=-1&amp;npep=0" title="An expression chip simulation of this model">chip</a> | 
<a href="/thegpm-cgi/phplctab.pl?path=$full_gpm_path&amp;proex=-1&amp;npep=0" title="A table of observed peptides">peptide</a> | 
<a href="/thegpm-cgi/ptable.pl?path=$full_gpm_path&sort=expect&amp;proex=-1&amp;npep=0" title="An exportable table of this model">table</a> | 
<a href="/thegpm-cgi/transform_xml.pl?type=all&file=..$full_gpm_path&amp;proex=-1&amp;npep=0" class="small_link" title="A detailed view of the model" target=_BLANK>details</a> | 
<a href="/thegpm-cgi/pgo.pl?path=$full_gpm_path&amp;proex=-1&amp;npep=0" class="small_link" title="GO classification diagram">GO</a> | 
<a href="/thegpm-cgi/ppath.pl?path=$full_gpm_path&amp;proex=-1&amp;npep=0" class="small_link" title="KEGG pathways analysis">path</a> | 
<a href="$full_gpm_path" title="The XML file that contains this model">XML</a> | 
<a href="http://$g_site{'wiki_url'}/wiki/$gpmnum" title="A wiki entry for the data set $gpmnum" target="_WIKI">wiki</a>
	);

	return $result;

} # end simple_model_toolbar

sub log_bad_gpm_file {
#
# This function will log the date/time of a request to what is deemed to be a corrupted gzip file along with the
#	path/filename.
#
# Inputs: the gpm filename and the path to the script that found the bad file
#
# Outputs: appends the current timestamp and GPM filename to the logfile specified below.
#

	my $LOG_NAME = '../gpm/corrupt_gzips.txt';
	my $gpmfile = shift;
	my $source_script = shift;

	if (not(-e $LOG_NAME)) { # create the log

		open(OUT, ">$LOG_NAME") or die("Error creating logfile: $!\n");

	} else { # open existing log

		open(OUT, ">>$LOG_NAME") or die("Error logging bad gzip: $!");

	} # end else
	
	print OUT localtime() . "\t$gpmfile\t$source_script\n";
	close(OUT);

} # end log_bad_gpm_file

sub unlikely {
	my ($_s,$_rx,$_h,$_l) = @_;
	my $s = $_s;
	$s =~ s/$_rx/$1,$2/g;
	%$_h = ();
	my @peptides = split /,/,$s;
	my $peptide;
	my $temp;
	my $pos = 0;
	my $length = 0;
	my $a;
	my $phobic;
	my $membrane;
	foreach $peptide(@peptides)	{
		$length = length($peptide);
		$temp = $peptide;
		$temp =~ s/[LIVMFWY]//g;
		$phobic = 1 - length($temp)/$length;
		$temp = $peptide;
		$temp =~ s/[LIVMACFWY]{7,}//g;
		if($peptide ne $temp)	{
			$membrane = 1;
		}
		else	{
			$membrane = 0;
		}

		if($length < 6)	{ 
			$a = 0;
			while($a < $length)	{
				$$_h{$a + $pos} = "$length 0.3";
				$a++;
			}
		}
		elsif($length > 35)	{ 
			$a = 0;
			while($a < $length)	{
				$$_h{$a + $pos} = "$length 0.4";
				$a++;
			}
		}
		elsif($length > 10 and $phobic >= 0.68)	{ 
			$a = 0;
			while($a < $length)	{
				$$_h{$a + $pos} = "$length 0.3";
				$a++;
			}
		}
		elsif($membrane)	{ 
			$a = 0;
			while($a < $length)	{
				$$_h{$a + $pos} = "$length 0.3";
				$a++;
			}
		}
		elsif(not ($_l =~ /gi\|/) and $peptide =~ /N[^P][ST]/)	{ 
			$a = 0;
			while($a < $length)	{
				$$_h{$a + $pos} = "$length 0.15";
				$a++;
			}
		}
		$pos += $length;
	}
}

sub unlikely_residues {
	my ($_s,$_rx,$_l) = @_;
	my $s = $_s;
	$s =~ s/$_rx/$1,$2/g;
	my @peptides = split /,/,$s;
	my $peptide;
	my $temp;
	my $length = 0;
	my $return = 0;
	my $phobic;
	my $membrane;
	foreach $peptide(@peptides)	{
		$length = length($peptide);
		$temp = $peptide;
		$temp =~ s/[LIVMFWY]//g;
		$phobic = 1 - length($temp)/$length;
		$temp = $peptide;
		$temp =~ s/[LIVMACFWY]{7,}//g;
		if($peptide ne $temp)	{
			$membrane = 1;
		}
		else	{
			$membrane = 0;
		}

		if($length < 6)	{ 
			$return += $length;
		}
		elsif($length > 35)	{ 
			$return += $length;
		}
		elsif($length > 10 and $phobic >= 0.68)	{ 
			$return += $length;
		}
		elsif(not ($_l =~ /gi\|/) and $peptide =~ /N[^P][ST]/)	{ 
			$return += $length;
		}
		elsif($membrane)	{ 
			$return += $length;
		}
	}
	return $return;
}

sub translate_label
{
	my ($_l,$_t) = @_;
	if($_t == 0)	{
		return $_l;
	}
	LoadHgnc($_l);
	if($_t == 1)	{
		if($g_ipi{$_l} =~ /^NP/)	{
			return $g_ipi{$_l};
		}
		else	{
			return $_l;
		}
	}
	elsif($_t == 2)	{
		if($g_hgnc{$_l} =~ /\w/)	{
			return $g_hgnc{$_l};
		}
		else	{
			return $_l;
		}
	}
	elsif($_t == 3)	{
		if($g_pubmed{$_l} =~ /\w/)	{
			return "ncbi:$g_pubmed{$_l}";
		}
		else	{
			return $_l;
		}
	}
	elsif($_t == 4)	{
		if($g_swiss{$_l} =~ /\w/)	{
			return $g_swiss{$_l};
		}
		else	{
			return $_l;
		}
	}
	return $_l;
}

sub common_normal_bar
{
	my ($npep,$url,$proex,$ltype,$path) = @_;
	my $b = qq(<a href="/thegpm-cgi/plist.pl?npep=$npep&path=$url&amp;proex=$proex&amp;ltype=$ltype" class="small_link" title="Summary navigation page for this model">model</a> | );
	$b .= qq(<a href="/thegpm-cgi/pintersect.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;ltype=$ltype" title="Locate all models containing a set of proteins">context</a> | );
	$b .= qq(<a href="/thegpm-cgi/pgroup.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;ltype=$ltype" title="Protein grouping based on peptide homology">group</a> | );
	$b .= qq(<a href="/thegpm-cgi/pgel.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;ltype=$ltype" title="A 1D and 2D PAGE simulation of this model">gel</a> | );
	$b .= qq(<a href="/thegpm-cgi/pchip.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;ltype=$ltype" title="An expression chip simulation of this model\">chip</a> | );
	$b .= qq(<a href="/thegpm-cgi/phplctab.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;ltype=$ltype" title="A table of observed peptides">peptide</a> | );
	$b .= qq(<a href="/thegpm-cgi/ptable.pl?path=$url&amp;proex=$proex&amp;npep=$npep&amp;ltype=$ltype" title="A table of observed proteins">table</a> | );
	$b .= qq(<a href="/thegpm-cgi/transform_xml.pl?type=all&file=$path&amp;proex=$proex&amp;npep=$npep&amp;ltype=$ltype" class="small_link" title="A detailed view of the model" target=_BLANK>details</a> | );
	return $b;
}

sub decorate
{
	my ($line) = @_;
	$line =~ s/\n/<br>/g;
	if($line  =~ /tranche key/i)	{
		my ($key) = $line  =~ /tranche key[\:\=]*\s*(\S+)\s*/im;
		$key =~ s/\<br\>//gm;
		$key = uri_escape($key);
		$line  =~ s/tranche key[\:\=]*\s*.(\S+)(\s*)/<a href\=\"https:\/\/proteomecommons\.org\/dataset\.jsp\?i=$key\" target=\"_TRANCHE\">Tranche&nbsp;<img src=\"\/pics\/tranche.png\" border=\"0\" \/><\/a>$2/i;
	}
	if($line  =~ /trache key/i)	{
		my ($key) = $line  =~ /trache key[\:\=]*\s*(\S+)\s*/im;
		$key =~ s/\<br\>//gm;
		$key = uri_escape($key);
		$line  =~ s/trache key[\:\=]*\s*.(\S+)(\s*)/<a href\=\"https:\/\/proteomecommons.org\/dataset\.jsp\?i=$key\" target=\"_TRANCHE\">Tranche&nbsp;<img src=\"\/pics\/tranche.png\" border=\"0\" \/><\/a>$2/i;
	}
	if($line  =~ /RANCHE KEY/i)	{
		my ($key) = $line  =~ /RANCHE KEY[\:\=]*\s*(\S+)\s*/im;
		$key =~ s/\<br\>//gm;
		$key = uri_escape($key);
		$line  =~ s/RANCHE KEY[\:\=]*\s*.(\S+)(\s*)/<a href\=\"https:\/\/proteomecommons.org\/dataset\.jsp\?i=$key\" target=\"_TRANCHE\">Tranche&nbsp;<img src=\"\/pics\/tranche.png\" border=\"0\" \/><\/a>$2/i;
	}
	if($line  =~ /TRNACHE KEY/i)	{
		my ($key) = $line  =~ /TRNACHE KEY[\:\=]*\s*(\S+)\s*/im;
		$key =~ s/\<br\>//gm;
		$key = uri_escape($key);
		$line  =~ s/TRNACHE KEY[\:\=]*\s*.(\S+)(\s*)/<a href\=\"https:\/\/proteomecommons.org\/dataset\.jsp\?i=$key\" target=\"_TRANCHE\">Tranche&nbsp;<img src=\"\/pics\/tranche.png\" border=\"0\" \/><\/a>$2/i;
	}
	if($line  =~ /Tranche data set 1660/)	{
		my $key = "H3tDepQtchkesY4nX/GpgZfKTwoCaOIharp/Tgw/CJBXTR2MMNI103n/nv8PGy8VuC4VQ5CiUwFIBh/n48ITOo+iMAMAAAAAAARaDw==";
		$key = uri_escape($key);
		my $value = qq(<a href="https://proteomecommons.org/dataset.jsp?i=$key" target="_TRANCHE">Tranche&nbsp;<img src="/llama-magic-html/tranche.png" border="0" /></a>, published in Global survey of human T leukemic cells by integrating proteomics and transcriptomics profiling, Mol Cell Proteomics. 2007 6:1343-53 <a href="http://www.ncbi.nlm.nih.gov/pubmed/17519225" target="_TRANCHE">PubMed</a>);
		$line  =~ s/Tranche data set 1660/$value $2/i;
	}
	if($line =~ /PRIDE ID[\:\=\s]/i)	{
		my ($key) = $line  =~ /pride id[\:\=\s]*\s*(\d+)/im;
		if($key)	{
			$line  =~ s/pride id[\:\=\s]*\s*(\d+)/<a href\=\"http:\/\/www.ebi.ac.uk\/pride\/directLink.do?experimentAccessionNumber=$key\" target=\"_TRANCHE\">PRIDE $key&nbsp;<img src=\"\/pics\/pride_logo.png\" border=\"0\" \/><\/a>/gmi;
		}
	}
	if($line =~ /PubMed ID[\:\=]\s*\d+/i)	{
		my ($key) = $line  =~ /PubMed ID[\:\=]\s*(\d+)/im;
		if($key)	{
			$line  =~ s/PubMed ID[\:\=]\s*(\d+)/\(<a href\=\"http:\/\/www\.ncbi\.nih\.gov\/pubmed\/$1\" target=\"_TRANCHE\">PubMed<\/a>\)/gmi;
		}
	}
	if($line =~ /pae\d{6}/i)	{
		$line  =~ s/(pae\d{6})/<a href\=\"http:\/\/www.peptideatlas.org\/repository\/\" target=\"_TRANCHE\">$1 <img src\=\"\/pics\/isb.gif\" border\=\"0\" \/><\/a>/gi;
	}
	if($line =~ /gpm\d{11}/i)	{
		$line  =~ s/(gpm\d{11})/<a href\=\"http:\/\/db2.thegpm.org\/thegpm\-cgi\/dblist_gpmnum.pl?gpmnum=$1\" target=\"_TRANCHE\">$1<\/a>/gi;
	}
	if($line =~ /PRIDE KEY[\:\=]/i)	{
		my ($key) = $line  =~ /pride key[\:\=]*\s*(\d+)/im;
		if($key)	{
			$line  =~ s/pride id[\:\=]*\s*(\d+)/<a href\=\"http:\/\/www.ebi.ac.uk\/pride\/directLink.do?experimentAccessionNumber=$key\" target=\"_TRANCHE\">PRIDE $key&nbsp;<img src=\"\/pics\/pride_logo.png\" border=\"0\" \/><\/a>/i;
		}
	}
	elsif($line =~ /PRIDE accession/i)	{
		my ($key) = $line  =~ /PRIDE accession\s*(\d+)/im;
		if($key)	{
			$line  =~ s/PRIDE accession\s*(\d+)/<a href\=\"http:\/\/www.ebi.ac.uk\/pride\/directLink.do?experimentAccessionNumber=$key\" target=\"_TRANCHE\">PRIDE $key&nbsp;<img src=\"\/pics\/pride_logo.png\" border=\"0\" \/><\/a>/i;
		}
	}
	if($line =~ /HuPA\_\d+/i)	{
		my ($key) = $line=~ /Hupa\_(\d+)/i;
		$line =~ s/HuPA\_\d+/<a href='http:\/\/www.humanproteinpedia.org\/data_display?exp_id=$key' target='HuPA'>HuPA\_$key<\/a>/i;
	}
	while($line =~ /\[\[\S+ .+\]\]/i)	{
		$line =~ s/\[\[(\S+) (.+?)\]\]/<a href=\"$1\" target=\"_OTHER\">$2<\/a>/g;
	}
	$line =~ s/(BTO\:\d+)/\<a href\=\"http\:\/\/www.ebi.ac.uk\/ontology\-lookup\/browse.do\?ontName\=BTO\&termId\=$1\" target\=\"_ONTO\"\>$1&nbsp;\<img src\=\"\/pics\/help.gif\" border=\"0\"\>\<\/a\>/g;
	$line =~ s/(MS\:\d+)/\<a href\=\"http\:\/\/www.ebi.ac.uk\/ontology\-lookup\/browse.do\?ontName\=PSI-MS\&termId\=$1\" target\=\"_ONTO\"\>$1&nbsp;\<img src\=\"\/pics\/help.gif\" border=\"0\"\>\<\/a\>/g;
	$line =~ s/(DOID\:\d+)/\<a href\=\"http\:\/\/www.ebi.ac.uk\/ontology\-lookup\/browse.do\?ontName\=DOID\&termId\=$1\" target\=\"_ONTO\"\>$1&nbsp;\<img src\=\"\/pics\/help.gif\" border=\"0\"\>\<\/a\>/g;
	$line =~ s/(CL\:\d+)/\<a href\=\"http\:\/\/www.ebi.ac.uk\/ontology\-lookup\/browse.do\?ontName\=CL\&termId\=$1\" target\=\"_ONTO\"\>$1&nbsp;\<img src\=\"\/pics\/help.gif\" border=\"0\"\>\<\/a\>/g;
	$line =~ s/(GO\:\d+)/\<a href\=\"http\:\/\/www.ebi.ac.uk\/ontology\-lookup\/browse.do\?ontName\=GO\&termId\=$1\" target\=\"_ONTO\"\>$1&nbsp;\<img src\=\"\/pics\/help.gif\" border=\"0\"\>\<\/a\>/g;
	$line =~ s/http\:\/\/www\.ncbi\.nlm\.nih\.gov\/peptidome\/search\/index\.shtml\?sample\=/http\:\/\/www\.ncbi\.nlm\.nih\.gov\/peptidome\/repository\//g;
	$line =~ s/http\:\/\/www\.ncbi\.nlm\.nih\.gov\/peptidome\/repository\/(PSM\d\d\d\d)\/*/http\:\/\/www\.thegpm\.org\/Peptidome\/PSM\/$1\.txt\.html/g;
	$line =~ s/http\:\/\/www\.ncbi\.nlm\.nih\.gov\/peptidome\/repository\/(PSE\d\d\d)\/*/http\:\/\/www\.thegpm\.org\/Peptidome\/PSE\/$1\.txt\.html/g;
	$line =~ s/([\w\.\-]+\@[\w\-]+\.[\w\.\-]+)/<a href=\"mailto:$1\">$1 <img src=\"\/pics\/email.png\" border=\"0\" \/><\/a>/g;
	$line =~ s/(\{\{[a-z\/]+\})\)/$1\}/g;
	$line =~ s/(\{\{[a-z\/]+\})([^\}])/$1\}$2/g;
	$line =~ s/\{\{/\</g;
	$line =~ s/\}\}/\>/g;
	$line =~ s/`/\"/g;
	$line =~ s/javascript//gi;
	$line =~ s/\son[mfbc]\w*\s*\=//gi;
	if($line =~ /\>PubMed\</i)	{
		$line =~ s/\>PubMed\</\>PubMed&nbsp<img src=\"\/pics\/ncbi_logo.png\" border=\"0\" \/>\</i;
	}
	return $line;		
}

sub GetSponsored
{
	my ($_i) = @_;
	if(not $g_site{"h_channel"} =~ /GPM\d\d\d/)	{
		return 0;
	}
	my $line = "|";
	$line .= " <a href=\"";
	$line .= "http://shop.bachem.com/ep6sf/search.ep?keyWords=$_i&categoryId=";
	$line .= "\" style=\"background-color:FFDDDD\" target=\"_sponsors\">Bachem</a> |";
	$line .= " <a href=\"";
	$line .= "http://www.picosearch.com/cgi-bin/ts.pl?index=437733&opt=ANY&query=$_i";
	$line .= "\" style=\"background-color:FFDDDD\" target=\"_sponsors\">Epitomics</a> |";
	$line .= " <a href=\"";
	$line .= "http://www.invitrogen.com/site/us/en/home/Global/invitrogen-search-results.html?searchTerm=$_i&searchRows=15&searchTypes=meta.collection%3Amrdb%3A&zeroCategories=meta.collection%3Acmgtbuyproduct%3A&noPagingRows=20";
	$line .= "\" style=\"background-color:FFDDDD\" target=\"_sponsors\">Invitrogen</a> |";
	$line .= " <a href=\"";
	$line .= "http://www.prosci-inc.com/shop/search.php?mode=search&query=$_i";
	$line .= "\" style=\"background-color:FFDDDD\" target=\"_sponsors\">ProSci</a> |";
	my $form = qq(<form target="_sponsors" name="thermo_form" style="display:none" method="post" action="http://www.bioreagents.com/search/searchResults.cfm" name="searchForm">
	<input name="searchText" type="text" value="$_i" /><input type="submit" value="Search" /></form>);
	$line .= $form;
	$line .= " <a ";
	$line .= " style=\"background-color:FFDDDD\" onClick=\"thermo_form.submit()\">Thermo</a> |";
	return $line;
}

sub GetPsimod
{
	my ($_k) = @_;
	my $v = sprintf("%.4f",$_k);
	if($g_trans{"empty"})	{
	%g_trans = 	(
				"-17.0265"	=>	"Ammonia-loss",
				"-17.0266"	=>	"Ammonia-loss",
				"-17.0160"	=>	"Ammonia-loss",
				"-17.0000"	=>	"Ammonia-loss",
				"-18.0106"	=>	"Dehydrated",
				"-18.0110"	=>	"Dehydrated",
				"-18.0000"	=>	"Dehydrated",
				"0.9848"	=>	"Deamidated",
				"0.9840"	=>	"Deamidated",
				"0.9970"	=>	"Label:+1 n",
				"1.0000"	=>	"Deamidated",
				"1.9940"	=>	"Label:+2 n",
				"1.9941"	=>	"Label:+2 n",
				"2.0042"	=>	"Label:+2 Da",
				"2.0043"	=>	"Label:+2 Da",
				"2.9911"	=>	"Label:+3 n",
				"3.0000"	=>	"Label:+3 Da",
				"3.0188"	=>	"Label:+3 Da",
				"3.0200"	=>	"Label:+3 Da",
				"3.9881"	=>	"Label:+4 Da",
				"4.0251"	=>	"Label:+4 Da",
				"4.0252"	=>	"Label:+4 Da",
				"4.0084"	=>	"Label:+4 Da",
				"4.0085"	=>	"Label:+4 Da",
				"6.0201"	=>	"Label:+6 Da",
				"6.0138"	=>	"Label:+6 Da",
				"6.0000"	=>	"Label:+6 Da",
				"6.0202"	=>	"Label:+6 Da",
				"7.0000"	=>	"Label:+7 Da",
				"8.0000"	=>	"Label:+8 Da",
				"8.0100"	=>	"Label:+8 Da",
				"8.0141"	=>	"Label:+8 Da",
				"8.0142"	=>	"Label:+8 Da",
				"9.0000"	=>	"Label:+9 Da",
				"10.0000"	=>	"Label:+10 Da",
				"10.0142"	=>	"Label:+10 Da",
				"10.0082"	=>	"Label:+10 Da",
				"10.0083"	=>	"Label:+10 Da",
				"10.0209"	=>	"Label:+10 Da",
				"10.0627"	=>	"Label:+10 Da",
				"10.0628"	=>	"Label:+10 Da",
				"12.0000"	=>	"Carbon",
				"14.0000"	=>	"Methyl",
				"14.0157"	=>	"Methyl",
				"14.0160"	=>	"Methyl",
				"14.0156"	=>	"Methyl",
				"15.9949"	=>	"Oxidation",
				"15.9990"	=>	"Oxidation",
				"16.0000"	=>	"Oxidation",
				"17.0345"	=>	"Methyl:+3 Da",
				"24.9952"	=>	"NTCB",
				"27.9949"	=>	"Formyl",
				"28.0313"	=>	"Dimethyl",
				"31.9898"	=>	"Dioxidation",
				"32.0000"	=>	"Dioxidation",
				"32.0564"	=>	"Dimethyl:+4 Da",
				"21.9819"	=>	"Cation:Na",
				"21.9820"	=>	"Cation:Na",
				"36.0757"	=>	"Dimethyl:+8 Da",
				"36.0455"	=>	"Dimethyl:+8 Da",
				"36.0456"	=>	"Dimethyl:+8 Da",
				"37.9558"	=>	"Cation:K",
				"37.9559"	=>	"Cation:K",
				"42.0105"	=>	"Acetyl",
				"42.0106"	=>	"Acetyl",
				"42.0000"	=>	"Acetyl",
				"42.0469"	=>	"Trimethyl",
				"42.0470"	=>	"Trimethyl",
				"43.0058"	=>	"Carbamyl",
				"43.0060"	=>	"Carbamyl",
				"43.9898"	=>	"Carboxyl",
				"44.0000"	=>	"Carboxyl",
				"45.0294"	=>	"Acetyl:2H(3)",
				"45.9877"	=>	"Methylthio",
				"47.9847"	=>	"Trioxidation",
				"47.9444"	=>	"Delta:S(-1)Se(1)",
				"56.0262"	=>	"Propionyl",
				"57.0220"	=>	"Carbamidomethyl",
				"57.0210"	=>	"Carbamidomethyl",
				"57.0215"	=>	"Carbamidomethyl",
				"57.0000"	=>	"Carbamidomethyl",
				"58.0185"	=>	"CAM+Label:+ 1 Da",
				"58.0055"	=>	"Carboxymethyl",
				"58.0054"	=>	"Carboxymethyl",
				"58.0000"	=>	"Carboxymethyl",
				"59.0363"	=>	"Propionyl:+3 Da",
				"68.0626"	=>	"Piperidine",
				"68.0374"	=>	"IMID",
				"68.0375"	=>	"IMID",
				"72.0625"	=>	"IMID:+2 Da",
				"72.0626"	=>	"IMID:+2 Da",
				"70.0400"	=>	"Crotonaldehyde",
				"70.0418"	=>	"Crotonaldehyde",
				"70.0419"	=>	"Crotonaldehyde",
				"71.0371"	=>	"Propionamide",
				"72.0211"	=>	"Carboxyethyl",
				"79.9568"	=>	"Sulfo",
				"79.9663"	=>	"Phospho",
				"80.0000"	=>	"Phospho",
				"86.0003"	=>	"Malonyl",
				"86.0004"	=>	"Malonyl",
				"86.0000"	=>	"Malonyl",
				"87.9982"	=>	"Thioacyl",
				"87.9983"	=>	"Thioacyl",
				"100.0000"	=>	"Succinyl",
				"100.0160"	=>	"Succinyl",
				"104.0290"	=>	"Succinyl:+4 Da",
				"104.0294"	=>	"Succinyl:+4 Da",
				"104.0295"	=>	"Succinyl:+4 Da",
				"104.0410"	=>	"Succinyl:+4 Da",
				"104.0411"	=>	"Succinyl:+4 Da",
				"104.0412"	=>	"Succinyl:+4 Da",
				"105.0000"	=>	"ICPL",
				"105.0200"	=>	"ICPL",
				"105.0214"	=>	"ICPL",
				"105.0210"	=>	"ICPL",
				"105.0215"	=>	"ICPL",
				"109.0000"	=>	"ICPL:2H(4)",
				"109.0500"	=>	"ICPL:2H(4)",
				"109.0465"	=>	"ICPL:2H(4)",
				"109.0470"	=>	"ICPL:2H(4)",
				"109.0466"	=>	"ICPL:2H(4)",
				"111.0000"	=>	"ICPL:13C(6)",
				"111.0415"	=>	"ICPL:13C(6)",
				"111.0420"	=>	"ICPL:13C(6)",
				"111.0416"	=>	"ICPL:13C(6)",
				"125.0477"	=>	"Nethylmaleimide",
				"125.0480"	=>	"Nethylmaleimide",
				"114.0429"	=>	"GlyGly",
				"114.0430"	=>	"GlyGly",
				"119.0040"	=>	"Cysteinyl",
				"119.0041"	=>	"Cysteinyl",
				"125.8970"	=>	"Iodo",
				"251.7930"	=>	"Diodo",
				"377.6900"	=>	"Triodo",
				"595.6130"	=>	"Tetraiodo",
				"125.8966"	=>	"Iodo",
				"251.7932"	=>	"Diodo",
				"251.7933"	=>	"Diodo",
				"377.6899"	=>	"Triodo",
				"595.6128"	=>	"Tetraiodo",
				"144.1000"	=>	"iTRAQ",
				"144.1100"	=>	"iTRAQ",
				"144.1020"	=>	"iTRAQ",
				"144.1021"	=>	"iTRAQ",
				"144.0000"	=>	"iTRAQ",
				"140.0949"	=>	"mTRAQ",
				"140.0950"	=>	"mTRAQ",
				"148.1092"	=>	"mTRAQ:+8Da",
				"148.1090"	=>	"mTRAQ:+8Da",
				"145.0197"	=>	"CAMthiopropanoyl",
				"145.0200"	=>	"CAMthiopropanoyl",
				"442.2250"	=>	"ICAT-D",
				"442.2249"	=>	"ICAT-D",
				"224.1525"	=>	"TMT",
				"225.1558"	=>	"TMT2plex",
				"227.1269"	=>	"ICAT-C",
				"227.1270"	=>	"ICAT-C",
				"229.0140"	=>	"PyridoxalPhosphate",
				"229.1629"	=>	"TMT6plex",
				"229.1630"	=>	"TMT6plex",
				"304.2050"	=>	"iTRAQ8plex",
				"304.2053"	=>	"iTRAQ8plex",
				"304.2054"	=>	"iTRAQ8plex",
				"299.1667"	=>	"cysTMT",
				"299.1670"	=>	"cysTMT",
				"304.1772"	=>	"cysTMT6plex",
				"304.1770"	=>	"cysTMT6plex",
				"383.2281"	=>	"LRGG",
				"383.2280"	=>	"LRGG",
				"383.2300"	=>	"LRGG",
				"383.2000"	=>	"LRGG",
				"396.0838"	=>	"Bacillithiol",
				"396.0840"	=>	"Bacillithiol",
				"396.0839"	=>	"Bacillithiol",
				"572.1812"	=>	"TMPP-Ac"
			);
	}
	my $c = $g_trans{$v};
	if($c)	{
		return $c;
	}
#	print "$v<hr>";
	return $_k;
}

sub load_chr
{
	my ($_t,$_m) = @_;
	if($_t =~ /ENSP\d+/)	{
		if($$_m{"human"})	{
			return 0;
		}
		if(not -e "../annotation/human_chr.txt")	{
			return 0;
		}
		open(IN,"<../annotation/human_chr.txt");
		my @v = <IN>;
		close(IN);
		my $l;
		my $acc;
		my $content;
		$$_m{"human"} = 1;
		foreach $l(@v)	{
			chomp($l);
			($acc) = $l =~ /(^ENSP\d+)/;
			if($acc)	{
				($content) = $l =~ /ENSP\d+(\t.+)/;
				$$_m{$acc} = $content;
			}
		}
		return 1;
	}
	elsif($_t =~ /ENSMUSP\d+/)	{
		if($$_m{"mouse"})	{
			return 0;
		}
		if(not -e "../annotation/mouse_chr.txt")	{
			return 0;
		}
		open(IN,"<../annotation/mouse_chr.txt");
		my @v = <IN>;
		close(IN);
		my $l;
		my $acc;
		my $content;
		$$_m{"mouse"} = 1;
		foreach $l(@v)	{
			chomp($l);
			($acc) = $l =~ /(^ENSMUSP\d+)/;
			if($acc)	{
				($content) = $l =~ /ENSMUSP\d+(\t.+)/;
				$$_m{$acc} = $content;
			}
		}
		return 1;
	}
	elsif($_t =~ /ENSRNOP\d+/)	{
		if($$_m{"rat"})	{
			return 0;
		}
		if(not -e "../annotation/rat_chr.txt")	{
			return 0;
		}
		open(IN,"<../annotation/rat_chr.txt");
		my @v = <IN>;
		close(IN);
		my $l;
		my $acc;
		my $content;
		$$_m{"rat"} = 1;
		foreach $l(@v)	{
			chomp($l);
			($acc) = $l =~ /(^ENSRNOP\d+)/;
			if($acc)	{
				($content) = $l =~ /ENSRNOP\d+(\t.+)/;
				$$_m{$acc} = $content;
			}
		}
		return 1;
	}
	return 0;
}

sub zeta
{
	my ($_p,$_z,$_t) = @_;
	my (@l) = $_p =~ /[RHK]/g;
	my (@p) = $_p =~ /P/g;
	my $cannon = scalar(@l) + 1;
	if($_z > $cannon and scalar(@p))	{
		$cannon++;
	}
	if($_t)	{
		return sprintf("%i/%i",$_z,$cannon);
	}
	return $_z/$cannon;
}

sub FormatInfo
{
	my ($e,$n) = @_;
	if(not $n)	{
		$n = 101;
	}
	$e =~ s/\<br\>/&nbsp;<a onClick=\"toggleBox(\'desc$n\')\">(show domains)<\/a><div id=\'desc$n\' style=\'display:none\'><br>/;
	$e .= "\n</div>\n";
	return $e;
}

sub get_page_cache
{
	my ($_c,$_cvs) = @_;
	my $path = $$_c{"path"};
	if(not open(IN,"<$path"))	{
		$$_c{"status"} = "read failed";
		$$_c{"ok"} = 0;
		return 0;
	}
	my @k;
	my $good = 1;
	my $l;
	my $c = 0;
	while(<IN>)	{
		chomp($_);
		if(/^#/)	{
			@k = split /#/,$_;
 			$_ = <IN>;
			chomp($_);
			my @lines;
			while($_ and /\\\\$/)	{
				s/\\\\$//;
				push(@lines,$_);
 				$_ = <IN>;
				chomp($_);

			}
			push(@lines,$_);
			my $ref = $$_cvs{@k[1]};
			if($ref)	{
				my @x;
				my $y;
				@$ref = ();
				foreach $y(@lines)	{
					@x = split /’\t‘/,$y;
					push(@$ref,@x);
				}
				$l = scalar(@$ref);
				if(@k[2] != $l)	{
					$good = 0;
					$$_c{"elements"}{@k[1]} = "@k[2]!=$l";
				}
				else	{
					$$_c{"elements"}{@k[1]} = "@k[2]=$l";
				}
				$c++;
			}
			else	{
				$$_c{"elements"}{@k[1]} = "bad array reference";
				$good = 0;
			}
		}
	}
	close(IN);
	if($c != scalar(keys(%$_cvs)))	{
		$c -= scalar(keys(%$_cvs));
		$$_c{"status"} = "bad array keys $c";
		unlink($path);
		$$_c{"ok"} = 0;
		return 0;
	}
	if($good)	{
		$$_c{"status"} = "good";
		$$_c{"ok"} = 1;
		return 1;
	}
	$$_c{"status"} = "bad";
	$$_c{"ok"} = 0;
	unlink($path);
	return 0;	
}

sub write_page_cache
{
	my ($_c,$_cvs) = @_;
	my $path = $$_c{"path"};
	if(not -e $$_c{"dir"})	{
		mkdir($$_c{"dir"});
	}
	if(not open(OUT,">$path"))	{
		$$_c{"status"} = "write failed";
		$$_c{"ok"} = 0;
		return 0;
	}
	my @keys = keys(%$_cvs);
	my $k;
	my $l;
	my $v;
	foreach $k(@keys)	{
		my $ref = $$_cvs{$k};
		$l = scalar(@$ref);
		print OUT "#$k#$l\n";
		$$_c{"elements"}{$k} = $l;
		my $i;
		my $len = 1;
		my $tl = 0;
		my $tl_max = 200000;
		foreach $i(@$ref)	{
			if($len < $l)	{
				$tl += length($i);
				if($tl > $tl_max)	{
					print OUT "$i\\\\\n";
					$tl = 0;
				}
				else	{
					print OUT "$i’\t‘";
				}
			}
			else	{
				print OUT "$i";
			}
			$len++;
		}
		print OUT "\n";
	}
	close(LOG);
	close(OUT);
	$$_c{"status"} = "good";
	$$_c{"ok"} = 1;
	return 1;
}

sub dump_cache
{
	my ($_c) = @_;
	my $line = "\n<!== page cache status report\n\n";
	$line .= "status = " . $$_c{"status"} . "\n";
	$line .= "path = " . $$_c{"path"} . "\n";
	$line .= "dir = " . $$_c{"dir"} . "\n";
	$line .= "ok = " . $$_c{"ok"} . "\n";
	my $ref = $$_c{"elements"};
	my @ks = keys(%$ref);
	my $k;
	foreach $k(@ks)	{
		$line .= "element $k: " . $$ref{$k} . "\n";
	}
	$line .= "\n==>\n";
	$$_c{"html_status"} = $line;
}

sub GetProteinDescription
{
	my ($id,$_max) = @_;
	if(not $_max)	{
		$_max = 60;
	}
	my $v = GetInfoNew($id);
	$v =~ s/[\r\n].+//si;
	$v =~ s/\&.+?;/ /sg;
	$v =~ s/\<.*?\>//sgi;
	$v =~ s/Annotated domains.*$//sgi;
	$v =~ s/\"/—/sgi;
	$v =~ s/\'/—/sgi;
	$v =~ s/\[.*?\]//sgi;
	$v =~ s/\&/ /g;
	if(length($v) > $_max)	{
		$v =~ s/^(.{$_max}\S+).+/$1 .../;
	}
	return $v;
}	

sub GetInfoNew
{
	my ($_l) = @_;
	my $v = GetInfo($_l);
	$v =~ s/^\s+//;
	$v =~ s/\s+/ /g;
	LoadHgnc($_l);
	my $d = $g_hgnc{$_l};
	if($d and not ($v =~ /\s$d[\s\,]/i or $v =~ /^$d[\s\,]/i))	{
		if(length($d) > 1)	{
			$v = qq(<a href="http://gpmdb.thegpm.org/$d" target="_GPMDB">$d</a>, $v);
		}
		else	{
			$v = qq($d, $v);
		}
	}
	elsif($v =~ /^\S\S+\,/)	{
		$v =~ s/^(\S\S+)\,/\<a href\=\"http\:\/\/gpmdb\.thegpm\.org\/$1\" target\=\"\_GPMDB\"\>$1\<\/a\>,/;

	}
	if(not $v)	{
		$v = "no protein information available";
	}
	$v =~ s/^(.{0,20})hypothetical protein/$1\<del title\=\'(sic)\'>hypothetical\<\/del\> \<ins\>observed\<\/ins\> protein/gi;
	return $v;
}

sub mrm_model
{
	my ($_s) = @_;
	my $score = 1;
	my @v = split //,$_s;
	my %aa;
	my $s;
	foreach $s(@v)	{
		$aa{$s} = $aa{$s} + 1;
	}
	if($aa{'K'} + $aa{'H'} + $aa{'R'} > 1)	{
		return -1;
	}
	if(@v[0] eq 'Q' or @v[0] eq 'C' or @v[0] eq 'E' or @v[0] eq 'K' or @v[0] eq 'R' or @v[0] eq 'P')	{
		return -2;
	}
	if(length($_s) < 9 or length($_s) > 30)	{
		return -3;
	}
	if(not $_s =~ /[KR]$/)	{
		return -4;
	}
	if($_s =~ /NG/)	{
		return -5;
	}
	return $score;
}

sub markup_annotation
{
	my ($_a) = @_;
	$_a =~ s/\s+//g;
	my @mods = split /\,/,$_a;
	my $isG = 0;
	if(scalar(@mods))	{
		my $y;
		foreach $y(@mods)	{
			if($y =~ /\-1\@B/)	{
				$isG = 1;
				$_a =~ s/\-1\@B\,*//;
				next;
			}
			$y =~ s/\@.+//;
			my $x = GetPsimod($y);
			$y =~ s/\+/\\\+/g;
			$y =~ s/\-/\\\-/g;
			$y =~ s/\./\\\./g;
			$y =~ s/\[/\\\[/g;
			$y =~ s/\]/\\\]/g;
			$_a =~ s/$y\@/$x\@/;
		}
		$_a =~ s/\,$//;
		$_a =~ s/\,/\, /g;
	}
	if($isG)	{
		$_a .= " (g)";
	}
	return $_a;
}
