#!/usr/bin/perl

#located in offline/framework/frog/

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;

sub creatfilelists;
sub printtags;
sub printruns;

my $buildtag;
my $cdbtag;
my $version;
my $verbose;
my $runnumber;
my $runlist;
my $printruns;
my $printtags;

my $noprint;

GetOptions('build:s' => \$buildtag, 'cdb:s' => \$cdbtag, 'list:s' => \$runlist, 'printruns' => \$printruns, 'printtags' => \$printtags, 'run:i' => \$runnumber, "verbose" =>\$verbose, 'version:s' => \$version);

my %ignore_datasets = (
    "mdc2" => 1,
    "rawdata" => 1
    );

my $dbh = DBI->connect("dbi:ODBC:FileCatalog","argouser") || die $DBI::error;
$dbh->{LongReadLen}=2000; # full file paths need to fit in here

if (defined $printtags)
{
    printtags();
}

if (defined $printruns)
{
    printruns();
}
if ($#ARGV < 0)
{
    if (! defined $noprint)
    {
	print "usage:\n";
	print "CreateDstList.pl <dsttype1> <dsttype2> ...\n";
	print "parameters:\n";
	print "--build <build tag> (mandatory)\n";
	print "--cdb <cdb tag> (mandatory)\n";
	print "--list <file with list of runs>\n";
	print "--printruns: print existing runs (for piping into a runlist)\n";
	print "--printtags: print existing tags\n";
	print "--run <run number>\n";
	print "--version <version> (only if we have multiple versions for same build/cdb tag)\n";
	print "--verbose print stuff\n";
    }
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
if (! defined $runnumber && ! defined $runlist)
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

my $getfiles =  $dbh->prepare("select filename from datasets where runnumber=? and dsttype=? and dataset = '$dataset' order by segment\n");
if (defined $runlist)
{
    if (! -f $runlist)
    {
	print "$runlist does not exist\n";
	exit(1);
    }
    open(F2,"$runlist");
    {
	while (my $run = <F2>)
	{
	    creatfilelists($run);
	}
    }
}
else
{
    creatfilelists($runnumber);
}

sub creatfilelists
{
    my $run = shift;
foreach my $ds (sort keys %dsttype)
{
    $getfiles->execute($run,$ds);
    if ($getfiles->rows == 0)
    {
	print "no run $run for dst type $ds and dataset $dataset\n";
	exit(1);
    }
    my $filename = sprintf("%s-%08d.list",lc $ds, $run);
    print "creating list for run $run --> $filename\n";
    open(F,">$filename");
    while (my @res = $getfiles->fetchrow_array())
    {
	print F "$res[0]\n";
    }
    close(F);
}
}

sub printtags
{
    my %dataset_exists = ();
    my $getdatasettypes = $dbh->prepare("select distinct(dataset) from datasets");
    $getdatasettypes->execute();
    while (my @res = $getdatasettypes->fetchrow_array())
    {
	if (! exists $ignore_datasets{$res[0]})
	{
	    $dataset_exists{$res[0]} = 1;
	}
    }
    foreach my $ds (sort keys %dataset_exists)
    {
	my @sp1 = split(/_/,$ds);
	print "build tag: $sp1[0], cdb tag: $sp1[1]";
	if ($#sp1 > 2)
	{
	    print "version: $sp1[2]";
	}
	print "\n";
    }
    exit(0);
}

sub printruns
{
    my $dataset;
    if (defined $buildtag && defined $cdbtag)
    {
	$buildtag =~ s/\.//g;
	$dataset = sprintf("%s_%s",$buildtag,$cdbtag);
	if (defined $version)
	{
	    $dataset = sprintf("%s_%s",$dataset,$version);
	}
    }
    if (defined $dataset)
    {
	my $getruns =  $dbh->prepare("select distinct(runnumber) from datasets where dataset = '$dataset' order by runnumber\n");
	$getruns->execute();
	if ($getruns->rows == 0)
	{
	    print "no run found for dataset $dataset\n";
	}
	else
	{
	    while (my @res = $getruns->fetchrow_array())
	    {
		print "$res[0]\n";
	    }
	}
    }
    else
    {
	my $sqlcmd = sprintf("select distinct(runnumber) from datasets where ");
	my $icnt = 0;
	foreach my $ignore (sort keys %ignore_datasets)
	{
	    if ($icnt > 0)
	    {
		$sqlcmd = sprintf("%s and ",$sqlcmd);
	    }
	    $sqlcmd = sprintf("%s dataset <> \'\%s\'",$sqlcmd,$ignore);
	    $icnt++;
	}
	$sqlcmd = sprintf("%s order by runnumber",$sqlcmd);
	my $getruns =  $dbh->prepare($sqlcmd);
	$getruns->execute();
	if ($getruns->rows == 0)
	{
	    print "no runs found\n";
	}
	else
	{
	    while (my @res = $getruns->fetchrow_array())
	    {
		print "$res[0]\n";
	    }
	}

    }
    exit(0);
}
