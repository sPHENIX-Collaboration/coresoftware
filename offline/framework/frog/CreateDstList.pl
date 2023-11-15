#!/usr/bin/perl

#located in offline/framework/frog/

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;

my $buildtag;
my $cdbtag;
my $version;
my $verbose;
my $runnumber;

GetOptions('build:s' => \$buildtag, 'cdb:s' => \$cdbtag, 'run:i' => \$runnumber, 'version:s' => \$version, "verbose" =>\$verbose);

if ($#ARGV < 0)
{
    print "usage:\n";
    print "CreateDstList.pl <dsttype1> <dsttype2> ...\n";
    print "parameters:\n";
    print "--build <build tag> (mandatory)\n";
    print "--cdb <cdb tag> (mandatory)\n";
    print "--run <run number> (mandatory)\n";
    print "--version <version> (only if we have multiple version for same build/cdb tag\n";
    print "--verbose print stuff\n";
    exit(0);
}
my %dsttype = ();
while($#ARGV >= 0)
{
  $dsttype{$ARGV[0]} = 1;
  shift (@ARGV);
}

if (! defined $buildtag)
{
    print "set build tag with --build <tag>\n";
    exit(1);
}
if (! defined $cdbtag)
{
    print "set cdb tag with --cdb <tag>\n";
    exit(1);
}
if (! defined $runnumber)
{
 print "set run number with --run <run number>\n";
 exit(1);
}
 

$buildtag =~ s/\.//g;

my $dataset = sprintf("%s_%s",$buildtag,$cdbtag);
if (defined $version)
{
    $dataset = sprintf("%s_%s",$dataset,$version);
}

my $dbh = DBI->connect("dbi:ODBC:FileCatalog","argouser") || die $DBI::error;
$dbh->{LongReadLen}=2000; # full file paths need to fit in here

my %dataset_exists = ();
my $getdatasettypes = $dbh->prepare("select distinct(dataset) from datasets");
$getdatasettypes->execute();
while (my @res = $getdatasettypes->fetchrow_array())
{
    $dataset_exists{$res[0]} = 1;
}
if (! exists $dataset_exists{$dataset})
{
    print "dataset $dataset does not exist, here is a list of our datasets\n";
    print "datasets are named buildtag_cdbtag_version\n";
    foreach my $ds (sort keys %dataset_exists)
    {
	print "$ds\n";
    }
    exit(1);
}
#$getdatasettypes->close();
my %dsttype_exists = ();
my $getdsttypes = $dbh->prepare("select distinct(dsttype) from datasets where dataset = '$dataset'");
$getdsttypes->execute();
while (my @res = $getdsttypes->fetchrow_array())
{
    $dsttype_exists{$res[0]} = 1;
}
foreach my $ds (sort keys %dsttype)
{
    if (! exists $dsttype_exists{$ds})
    {
	print "dst type $ds does not exist, here is a list of our dst types for dataset $dataset\n";
	foreach my $dstt (sort keys %dsttype_exists)
	{
	    print "$dstt\n";
	}
	exit(1);
    }
}

my $getfiles =  $dbh->prepare("select filename from datasets where runnumber=$runnumber and dsttype=? and dataset = '$dataset' order by segment\n");
foreach my $ds (sort keys %dsttype)
{
    $getfiles->execute($ds);
    if ($getfiles->rows == 0)
    {
	print "no run $runnumber for dst type $ds and dataset $dataset\n";
	exit(1);
    }
    my $filename = sprintf("%s-%08d.list",lc $ds, $runnumber);
    open(F,">$filename");
    while (my @res = $getfiles->fetchrow_array())
    {
	print F "$res[0]\n";
    }
    close(F);
}
