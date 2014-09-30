#!C:/Perl/bin/perl 

#    llama_magic_score.pl
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

#web interface for programs to find the best Llama nanobodies based on peptide mapping (from MS data) and CDR regions
#INPUT: DNA sequence files and MS files, XTandem parameters
#OUTPUT: HTML page with ranked, grouped nanobodies (also showing peptide mapping and highlighted CDR1/2/3 regions)
#DATA FLOW: DNA seqeunces -> Protein Sequences (longest non-reduncdant reding frame is used) -> digested with Trypsin ->
#XTandem search with MS data against the pre-digested protein sequences -> Map XTandem peptides (> min. Expectation)
#back to the sequences -> find CDR regions -> rank nanobodies based on coverage of CDR regions, etc. -> group nanobodies
#with the same CDR regions (0 or 1 difference in AA) -> HTML displayed output of this list

use strict;
use warnings;
use CGI ':standard';
use Proc::Background;
#use File::Path qw(make_path remove_tree);
#changed to not use this package since old version of ActiveState Perl on Rockefeller server and ppm not supported to get packages any longer 

my $DBLIST_FILE = "db_list.txt";
my $MSLIST_FILE = "ms_list.txt";
my $DBLIST_VIEWFILE = "db_list.txt";
my $MSLIST_VIEWFILE = "ms_list.txt";
my $SHOW_SCORE = 1;
my $SETTINGS_FILE = "../settings.txt"; 
my $img_source_plus = '/llama-magic-html/plus.gif';
my $img_source_minus = '/llama-magic-html/minus.gif';
my $img_source_star = '/llama-magic-html/greyplus.gif';
my $img_source_new = '/llama-magic-html/add_item.png';

my $DEVELOPER_VERSION = 1;
#my $DEVELOPER_VERSION = 0;
my $DEVELOPER_LOGFILE = "cgi_dev_log.txt";
my $BASE_DIR = "C:/NCDIR/Llama"; #default, changed when settings.txt file is read 
my $RESULTS_DIR = "results"; 
my %ALLOWED_DNA_FILE_TYPES = ('fa' => '1', 'fas' => '1', 'fast' => '1', 'fasta' => '1', 'fastq' => '1', 'fq' => '1');
my %ALLOWED_MS_FILE_TYPES = ('mgf' => '1');

eval #for exception handling
{
    if($DEVELOPER_VERSION) { open(DEVEL_OUT, ">>$BASE_DIR/$DEVELOPER_LOGFILE"); }
    if($DEVELOPER_VERSION) { print DEVEL_OUT "Opened developer log file...\n"; }
    
    my $err = "";
    if($err = read_settings())
    {
	display_error_page("Cannot load settings file: $err");
	if($DEVELOPER_VERSION) { print DEVEL_OUT "Cannot load settings file: $err\n"; }
    }
    else
    {
	if($DEVELOPER_VERSION) { print DEVEL_OUT "Loaded settings file: BASE_DIR = $BASE_DIR\n"; }
	if(!param())
	{#no posted data
	    display_home_page();
	}
	else
	{#form data was posted: process the data
	    my @params = param();
	    my $action = param('submit');
	    
	    if($DEVELOPER_VERSION) { print DEVEL_OUT "$action\n"; }
	    
	    if ($action eq 'Home')
	    {
		display_home_page();
	    }
	    elsif($action eq 'Upload')
	    {
		#get db name from html form
		my $db_name = param('db_name');
		if ($db_name eq "")
		{
		    display_error_page("Error: Please provide a name for the new database.");
		}
		else
		{
		    my $err_str = "";
		    
		    #create new db: creates directory structure, uploads files, and adds to db list
		    my $new_dir = "";
		    $err_str = create_new_db($db_name, $new_dir);
		    if ($err_str) { die $err_str; }
		    
		    #run translation from dna to protein and digest with trypsin -> protein db files
		    #(run asynchronously so the CGI program can return...)
		    #`"run_llama_scripts.pl" "db_scripts" "$BASE_DIR/$RESULTS_DIR/$new_dir"`;
		    my $proc1 = Proc::Background->new('C:/perl/bin/perl.exe', 'run_llama_scripts.pl', 'db_scripts', "$BASE_DIR/$RESULTS_DIR/$new_dir");
		    
		    #return informative message to user...
		    display_message_page("Your new database has been uploaded.  You should see it listed on the home page.  \
					 The link for the new database will become active once the predigested protein files have been created."); 
		}
	    }
	    elsif($action eq 'Create Candidate List')
	    {
		#get ms search name from html form
		my $db_id = param('ms_db_id');
		my $ms_name = param('ms_name');
		my $parent_err = param('mass_error');
		my $frag_err = param('frag_mass_error');
		if ($ms_name eq "") { display_error_page("Error: Please provide a name for this DB Search/Candidate List."); }
		elsif($parent_err eq "" or $parent_err <= 0) { display_error_page("Error: Parent mass error out of range."); }
		elsif($frag_err eq "" or $frag_err <= 0) { display_error_page("Error: Fragment mass error out of range."); }
		else
		{
		    my $err_str = "";
		    
		    #first, create new dirs and upload ms files to the current project
		    #also update db list txt file
		    my $new_dir = "";
		    $err_str = create_new_search($db_id, $ms_name, $new_dir);
		    if ($err_str) { die $err_str; }
		    
		    #(following will run asynchronously so the CGI program can return...)
		    #do tandem search, then, run map_peptides_to_protein and create candidate list html file
		    my $proc1 = Proc::Background->new('perl.exe', 'run_llama_scripts.pl', 'search_and_map_scripts', "$BASE_DIR/$RESULTS_DIR/$db_id/$new_dir", "$BASE_DIR/$RESULTS_DIR/$db_id", "$parent_err", "$frag_err", "/$RESULTS_DIR/$db_id/$new_dir", $SHOW_SCORE);
		    
		    #return informative message to user...
		    display_message_page("Your new search has been submitted.  You should see it listed on the home page under the database you selected.  \
					 The link for the Candidate List will become active once the XTandem! search is complete and the peptide mapping program has ranked the results."); 
		}
	    }
	}
    }
};
if ($@)
{
	if($DEVELOPER_VERSION) { print DEVEL_OUT "Exception thrown: $@\n"; }
	display_error_page("$@"); 
}

if($DEVELOPER_VERSION) { close(DEVEL_OUT); }

###########################

sub display_home_page
{
    my $type = shift;
    my $id;
    if ($type == 2) { $id = shift; }
    
    top_header(2);
    
    print "<tr class='main_list'>";

    #get list of dbs existing in repository any ms files + 'nanobody candidate results' under it
    open(IN, "$BASE_DIR/$RESULTS_DIR/$DBLIST_VIEWFILE");
    my @db_lines = <IN>;
    close(IN);
    
    open(IN, "$BASE_DIR/$RESULTS_DIR/$MSLIST_VIEWFILE");
    my @ms_lines = <IN>;
    close(IN);
    
    print "<td>",
	qq!<br/><b><u>Uploaded Databases:</u></b> <img src="$img_source_new" onclick="switch_right('db_form',0,0)" style="cursor:hand;" /><br /><br />!;
    for(my $i = 0; $i <= $#db_lines; $i++)
    {
	chomp($db_lines[$i]);
	if($db_lines[$i] =~ /^(\d+)\t([^\t]+)$/)
	{
	    my $id = $1; my $name = $2; 
	    my @ms_list;
	    #check status, if db translated and digested, create link, else if errors, create msg, else if not done, create disabled 'link' and msg
	    my $status = check_status_file("$id/protein"); #change to submit form!
	    if ($status eq 'DONE')
	    {
		#get ms searches for this db
		for(my $j = 0; $j <= $#ms_lines; $j++)
		{
		    if ($ms_lines[$j] =~ /^$id\t(\d+)\t([^\t]+)$/)
		    {
			push(@ms_list, [$1, $2]);
		    }
		}
		
		if ($#ms_list >= 0)
		{
		    print qq!<img src="$img_source_plus" onclick="ec('text_$id', 'img_$id')" id="img_$id" style="cursor:hand;" alt="+" />!;
		}
		else
		{
		    print qq!<img src="$img_source_star" alt="*" />!;
		}
		print qq!<a href="#" onclick="switch_right('ms_form', '$id', '$name')">$name</a><br />!;
	    }
	    elsif($status eq '') { print qq!<img src="$img_source_star" alt="*" />!, u($name), " &nbsp;(in progress...)", br(); }
	    else { print qq!<img src="$img_source_star" alt="*" />!,  u($name), " &nbsp;<a href='#' onclick='alert(\"" . $status . "\")'>*Error!</a>", br(), br(); } 
	    
	    #print the ms sub projects
	    print qq!<br/><div id="text_$id" style="display:none">!; 
	    for(my $j = 0; $j <= $#ms_list; $j++)
	    {
		my $ms_id = $ms_list[$j][0]; my $ms_name = $ms_list[$j][1];
		
		print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
		
		#check status, if tandem run/candidate list created, create link, else if errors, create msg, else if not done, create disabled 'link' and msg
		my $status = check_status_file("$id/$ms_id");
		if ($status eq 'DONE')
		{
		    my @files = <$BASE_DIR/$RESULTS_DIR/$id/$ms_id/tandem/results/*.1.cdr_coverage.html>;
		    #strip everything but file name....
		    $files[0] =~ /[\/\\]([^\/\\]+)$/;
		    print a({href=>"../llama-magic/$id/$ms_id/tandem/results/$1", target=>'_blank'}, "$ms_name"), br(), br();
		}
		elsif($status eq '') { print u($ms_name), " &nbsp;(in progress...)", br(), br(); }
		else { print u($ms_name), " &nbsp;<a href='#' onclick='alert(\"" . $status . "\")'>*Error!</a>", br(), br(); } 
		
	    }
	    if($#ms_list == -1) { print br(); }
	    print '</div>';
	}
    }
    
    print "</td>";
    
    print "<td>";
    display_help();
    display_db_upload_form();
    display_ms_upload_form();
    print "</td>";
    
    print "</tr>";
    display_footer();
}

sub display_help
{
    print "<div id='help'>",
	p(b(u('Instructions:')));
    print <<HELPTEXT;
<p class="instruction_text">1. Upload the database (fasta) files.  (Use the plus sign to add a new database.)  <br /><br />
2. Upload the MS (mgf) files. (Click on the database name once it has been added to the list on the left.)<br><br>
3. The result will be the nanobody candidate list (html) file.  <br><br>
HELPTEXT
#<i>Please see the </i><b><u>About</u></b> <i>page for more information regarding the format of the result
#file and the overall procedure used. </i></p>

    print "</div>";
}

sub display_db_upload_form
{
    print qq!<hr /><div id="db_form" style="display:none;"><a id="close" href="#" onclick="close_right('db_form')"></a>!,
	p(b(u('Upload a new Sequence Database:'))),
	'Name: ',
	textfield(-name=>'db_name', -value=>'', -size=>20, -maxlength=>20),
	br(), br(),
	'Select files (fasta format):',
	br(), 
	'<input name="dna_files" type="file" class="multi"/>',
	br(),
	submit(-name=>'submit', -id=>'db_upload_button', -value=>'Upload'), 
	"</div>";
}

sub get_db_name
{
    my $id = shift;
    
    open(IN, "$BASE_DIR/$RESULTS_DIR/$DBLIST_FILE");
    my @lines = <IN>;
    close(IN);
    
    for(my $i = 0; $i <= $#lines; $i++)
    {
	chomp($lines[$i]);
	if ($lines[$i] =~ /^$id\t(.+)$/)
	{
	    return $1;
	}
    }
    
    return '';
}

sub display_ms_upload_form
{
    print qq!<hr /><div id="ms_form" style="display:none;"><a id="close" href="#" onclick="close_right('ms_form')"></a>!,
	p(b(u('Create a new candidate list from MS files:'))), 
	"DB Name: ",
	qq!<input type="text" name="ms_db_name" id="ms_db_name" value="" disabled)!, 
	br(), br(),
	'Name: ',
	textfield(-name=>'ms_name', -value=>'', -size=>20, -maxlength=>20),
	br(), br(), 
	'Select MS files (mgf format):',
	br(), 
	'<input name="ms_files" type="file" class="multi"/>',
	br(), 
	u('X! Tandem parameters:'), br(), 
	'Parent mass error: +/- ',
	textfield(-name=>'mass_error', -value=>'10', -size=>4, -maxlength=>4), ' ppm', br(),
	'Fragment mass error:&nbsp;',
	textfield(-name=>'frag_mass_error', -value=>'0.4', -size=>4, -maxlength=>4), ' Daltons', br(), br(), 
	hidden(-name=>'ms_db_id', -id=>'ms_db_id', -value=>""),
	submit(-name=>'submit', -id=>'ms_upload_button', -value=>'Create Candidate List'), 
	"</div>";
}

sub create_new_db
{##add feature to lock db_list.txt file?
    my $new_db_name = $_[0];
    
    #read in last db created from db list, and increment ID for dir name
    open(IN, "$BASE_DIR/$RESULTS_DIR/$DBLIST_FILE");
    my @lines = <IN>;
    close(IN);
    
    my $new_dir = 1;
    my $i;
    for($i = $#lines; $i >= 0; $i--)
    {
	if ($lines[$i] =~ /^(\d+)\t(.+)/){ last; }
    }
    if ($i >= 0)
    {
	$lines[$i] =~ /(\d+)\t(.+)/;
	$new_dir = $1+1;
    }
    
    #create new dir structure for the db
    my $err_str = "";
    
    #make_path("$BASE_DIR/$RESULTS_DIR/$new_dir/dna", "$BASE_DIR/$RESULTS_DIR/$new_dir/protein");
    mkdir("$BASE_DIR/$RESULTS_DIR/$new_dir");
    mkdir("$BASE_DIR/$RESULTS_DIR/$new_dir/dna");
    mkdir("$BASE_DIR/$RESULTS_DIR/$new_dir/protein");
    
    #upload dna fasta files:
    if (!$err_str) { $err_str = upload_files('dna', "$new_dir/dna"); }
    
    #create taxonomy file for X! Tandem
    if (!$err_str) { $err_str = create_taxonomy_xml("$new_dir"); }
    
    #if error during upload, remove newly created db 
    if ($err_str)
    {
	#remove_tree("$BASE_DIR/$RESULTS_DIR/$new_dir");
	unlink glob "$BASE_DIR/$RESULTS_DIR/$new_dir/*.*";
	unlink glob "$BASE_DIR/$RESULTS_DIR/$new_dir/dna/*.*";
	rmdir("$BASE_DIR/$RESULTS_DIR/$new_dir/dna");
	rmdir("$BASE_DIR/$RESULTS_DIR/$new_dir/protein");
	rmdir("$BASE_DIR/$RESULTS_DIR/$new_dir");
    }
    else
    {#no error, add db to db list txt file
	open(OUT, ">>$BASE_DIR/$RESULTS_DIR/$DBLIST_FILE");
	print OUT "$new_dir\t$new_db_name\n";
	close(OUT);
	
	#ALSO, if db list view file is different, add it there as well
	if (not($DBLIST_VIEWFILE eq $DBLIST_FILE))
	{
	    open(OUT, ">>$BASE_DIR/$RESULTS_DIR/$DBLIST_VIEWFILE");
	    print OUT "$new_dir\t$new_db_name\n";
	    close(OUT);
	}
	
	$_[1] = $new_dir;
    }
    return $err_str;
}

sub create_new_search
{##add feature to lock ms_list.txt file?
    my $db_id = $_[0];
    my $new_ms_name = $_[1];
    
    #read list of db searches, get next search id for the db id selected
    open(IN, "$BASE_DIR/$RESULTS_DIR/$MSLIST_FILE");
    my @lines = <IN>;
    close(IN);
    
    my $new_ms_id = 1;
    for(my $i = 0; $i <= $#lines; $i++)
    {
	if ($lines[$i] =~ /(\d+)\t(\d+)\t(.+)/)
	{
	    if ($1 eq $db_id and $2 >= $new_ms_id) { $new_ms_id = $2+1; }
	}
    }
    
    #create new dir structure for the search
    my $err_str = "";
    
    #make_path("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/mgf", "$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/tandem/results");
    mkdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id");
    mkdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/mgf");
    mkdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/tandem");
    mkdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/tandem/results");
    
    #upload mgf files:
    if (!$err_str) { $err_str = upload_files('ms', "$db_id/$new_ms_id/mgf") }
    
    if ($err_str)
    {
	#if error during upload, remove newly created dir and search files
	#remove_tree("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id");
	unlink glob "$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/*.*";
	unlink glob "$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/mgf/*.*";
	rmdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/tandem/results");
	rmdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/tandem");
	rmdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id/mgf");
	rmdir("$BASE_DIR/$RESULTS_DIR/$db_id/$new_ms_id");
    }
    else
    {#if no error, add the search to the search list file
	
	open(OUT, ">>$BASE_DIR/$RESULTS_DIR/$MSLIST_FILE");
	print OUT "$db_id\t$new_ms_id\t$new_ms_name\n"; 
	close(OUT);
	
	#ALSO, if MS list view file is different, add it there as well
	if (not($MSLIST_VIEWFILE eq $MSLIST_FILE))
	{
	    open(OUT, ">>$BASE_DIR/$RESULTS_DIR/$MSLIST_VIEWFILE");
	    print OUT "$db_id\t$new_ms_id\t$new_ms_name\n"; 
	    close(OUT);
	}
	
	$_[2] = $new_ms_id;
    }
    return $err_str;
}

sub display_footer
{
	print '</table>',
	      end_multipart_form(),
	      end_html();
}

sub display_error_page
{
    my $msg = shift;
    top_header(1);
    print "<tr class='main_help'><td>";
    print p($msg);
    print "</td></tr>";
    
    display_footer();
}

sub display_message_page
{
    my $msg = shift;
    top_header(1);
    print "<tr class='main_help'><td>";
    print p($msg),
    a({href=>"../llama-magic-cgi/llama_magic_score.pl?submit=Home"}, "Home");
    print "</td></tr>";
    display_footer();
}

#########################

sub read_settings
{
	open(IN, "$SETTINGS_FILE") || return $!;
	my $found = 0;
	while(<IN>)
	{
		chomp();
		if(/^INSTALL_DIR=(.*)$/) { $BASE_DIR = $1; $found++; }
	}
	close(IN);
	if($found == 1) { return ""; }
	else { return "Information missing from settings file: $SETTINGS_FILE.\n"; }
}

sub top_header
{
    my $n = shift;
    print header(),
	  start_html(-title => 'Llama Magic',
		     -style => [{ -src=>'../llama-magic-html/smoothness/jquery-ui-1.10.3.custom.min.css' }, #download and fix this
				{ -src=>'../llama-magic-html/main.css' }], 
		     -script => [{ -type=>'javascript', -src=>'../llama-magic-html/jquery-1.9.1.min.js' },
				{ -type=>'javascript', -src=>'../llama-magic-html/jquery-ui-1.10.3.custom.min.js' }, #download and fix this
				{ -type=>'javascript', -src=>'../llama-magic-html/jquery.MultiFile.pack.js' },
				{ -type=>'javascript', -src=>'../llama-magic-html/setup.js' }]);
	  print qq!<div id="dialog" title="Uploading..."><p>Your files are being uploaded to the server.  Please be patient, the page will refresh when the upload has completed...</p></div>!;
	  print start_multipart_form(-method=>'POST', -action=>"../llama-magic-cgi/llama_magic_score.pl"), 
	  "<table class='maintable'><tr class='banner' ><td colspan='$n'>",
	  qq!<h1><img id="llama" src="/llama-magic-html/llama.jpg"> Welcome to Llama Magic </h1>!,
	  #a({href=>"../llama-magic-cgi/llama_magic_score.pl?submit=Home"}, "Home"),
	  #  "&nbsp;|&nbsp;",
	  #  a({href=>"../llama-magic-html/About-Llama-Magic.htm", target=>'_blank'}, "About"),
	    "</td></tr>";
}

sub upload_file
{
	my $user_fname = $_[0];
	my $fh = $_[1];
	my $local_fname_root = $_[2];

	my $extension = $user_fname;
	if ($extension =~ s/^.*\.([^\.]+)$/$1/)
	{
		my $line;
		if(open (OUTFILE, ">$local_fname_root.$extension"))
		{
			#save the file, create the local name using the extension of the user file
			if ($user_fname =~ /\.txt$/i)
			{
				while ( $line=<$fh> )
				{
					chomp($line);
					$line=~s/\r$//;
					$line=~s/\r([^\n])/\n$1/g;
					print OUTFILE "$line\n";
				}
			}
			else
			{
				binmode OUTFILE;
				while ($line=<$fh>)
				{
					print OUTFILE $line;
				}
			}
			close(OUTFILE);
		}
		else { return  "Could not open local file to save $user_fname.\n"; }
	}
	else { return "Could not extract extension from $user_fname.\n"; }

	$_[3] = $extension;
	return "";
}

sub upload_files
{
    my $file_type = shift;
    my %allowed_exts;
    if ($file_type eq 'dna') { %allowed_exts = %ALLOWED_DNA_FILE_TYPES; }
    elsif($file_type eq 'ms') { %allowed_exts = %ALLOWED_MS_FILE_TYPES; }
    
    my $new_dir = shift;
    my $err_str = "";
    
    #upload files:
    my @lw_fh = upload("$file_type" . '_files'); # undef may be returned if it's not a valid file handle, e.g. file transfer interrupted by user
    if (!@lw_fh || $#lw_fh == -1)
    {
	return "Error: No files uploaded.";
    }
    
    my @remote_files = param("$file_type" . '_files');
    my $lw_fh; my $i = 1;
    foreach $lw_fh (@lw_fh)
    {
	    if (defined $lw_fh)
	    {
		    my $local_fname = "$BASE_DIR/$RESULTS_DIR/$new_dir/$file_type" . "_$i";
		    my $io_fh = $lw_fh -> handle; # Upgrade the handle to one compatible with IO::Handle:
		    my $ext;
		    if($err_str = upload_file($remote_files[$i-1], $io_fh, $local_fname, $ext))
		    { last; }
		    
		    if(!defined $allowed_exts{lc $ext})
		    {
			    my $ext_list = join ', ', keys %allowed_exts;
			    $err_str = "Unrecognized file type: $remote_files[$i-1].\n";
			    last;
		    }
		    my $remote_fname_root = $remote_files[$i-1];
		    $remote_fname_root =~ s/\.\w\w\w$//; #remove extension
		    $i++;
	    }
	    else { $err_str = "Could not upload file.\n"; last; }
    }
    return $err_str;
}

sub create_taxonomy_xml
{
    my $dir = shift;
    if(!open(OUT, ">$BASE_DIR/$RESULTS_DIR/$dir/taxonomy.xml")) { return "Error creating taxonomy file: '$BASE_DIR/$RESULTS_DIR/$dir/taxonomy.xml'"; }
    print OUT qq(<?xml version="1.0"?><bioml label="x! taxon-to-file matching list">\
		<taxon label="llama"><file format="peptide" URL="../results/$dir/protein/longest_nr_predigested.fasta" />
		</taxon></bioml>);
    close(OUT);
    return "";
}

sub check_status_file
{
    my $dir = shift;

    #try to open file:
    if (-e "$BASE_DIR/$RESULTS_DIR/$dir/status.txt")
    {
	open(IN, "$BASE_DIR/$RESULTS_DIR/$dir/status.txt");
	my @lines = <IN>;
	close(IN);
	my $err_str = '';
	my $done = 0;
	for(my $i = 0; $i <= $#lines; $i++)
	{
	    chomp($lines[$i]);
	    if ($lines[$i] =~ /^ERROR:/)
	    {
		if ($lines[$i] =~ /exited with value/) 
		{
		    chomp($lines[$i+1]);
		    if($lines[$i+1] ne ''){ $err_str .= $lines[$i] . ': ' . $lines[$i+1]; }
		    else { $err_str .= $lines[$i] }
		}
		else { $err_str .= $lines[$i]; }
		$err_str .= '\n';
	    }
	    elsif ($lines[$i] eq 'DONE')
	    {
		$done = 1;
		last;
	    }
	}
	if ($done)
	{
	    if ($err_str) { return $err_str; }
	    else { return 'DONE'; }
	}
	else { return ""; }
    }
    else { return ""; }
}