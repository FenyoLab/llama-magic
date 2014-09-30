#!/usr/local/bin/perl

## Version 2004.06.11
## Version 2004.08.06
## Version 2004.10.18 - added get_root() and get_cache_root()
## defines.pl
## Copyright (C) 2004 Ronald C Beavis, all rights reserved
## The Global Proteome Machine 
## This software is a component of the X! proteomics software
## development project
##
## Use of this software governed by the Artistic license,
## as reproduced at http://www.opensource.org/licenses/artistic-license.php
##

return 1;
use strict;

sub get_server_name
{
	return "10.193.36.219:8082";
}

sub get_gpm_number
{
	return "GPM222";
}

sub get_root
{
	return "..";
}

sub get_cache_root
{
	return "..";
}

sub get_proxy
{
# Return an empty string if there is no proxy server

	return "";

# To include a proxy server, comment out the line above and insert
# the ip address or domain name of your proxy server in the line 
# below, and replace "192.168.0.1" with that value. Most proxy servers
# use port 8080, but this may be different in some environments.

	return "http://192.168.0.1:8080";
}
