#!/usr/bin/env perl

#located in offline/framework/frog/

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;

sub creatfilelists;
sub printtags;
sub printruns;

my $tag;
my $verbose;
my $runnumber;
my $runlist;
my $printruns;
my $printtags;
my $hpss;
my $noprint;

GetOptions('hpss' =>\$hpss, 'list:s' => \$runlist, 'printruns' => \$printruns, 'printtags' => \$printtags, 'run:i' => \$runnumber, 'tag:s' => \$tag, "verbose" =>\$verbose);

my %ignore_datasets = (
    "mdc2" => 1,
    );

my $dbh = DBI->connect("dbi:ODBC:FileCatalog_read") || die $DBI::error;
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
	print "--tag <tag (build_cdb_version)> (mandatory)\n";
	print "--hpss <print files in hpss>\n";
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

if (! defined $tag)
{
    print "set tag with --tag <tag>\n";
    exit(1);
}
if (! defined $runnumber && ! defined $runlist)
{
    print "set run number with --run <run number>\n";
    exit(1);
}

my $sqlcmd = sprintf("select filename,dsttype,segment from datasets where runnumber=? and tag = '$tag'");
if (defined $hpss)
{
    $sqlcmd = sprintf("%s and status = 0",$sqlcmd)
}
else
{
    $sqlcmd = sprintf("%s and status > 0",$sqlcmd)
}
my $getfiles =  $dbh->prepare($sqlcmd);
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
    my %files = ();
    my $run = shift;
    $getfiles->execute($run);
    while(my @res = $getfiles->fetchrow_array)
    {
	$files{$res[2]}{$res[1]}{$res[0]} = $res[1];
    }
    foreach my $dsttyp (keys %dsttype)
    {
	my $filename = sprintf("%s-%08d.list",lc $dsttyp, $run);
	open(F,">$filename");
	my $fcnt = 0;
	foreach my $segment (sort { $a <=> $b } keys %files)
	{
	    foreach my $file (keys %{$files{$segment}{$dsttyp}})
	    {
		print F "$file\n";
		$fcnt++;
	    }
	}
	close(F);
	if ($fcnt == 0)
	{
	    print "no dsts for type $dsttyp for run $run\n";
	    unlink $filename;
	}
    }
    #x	print Dumper(\%files);
}

sub printtags
{
    my $sqlcmd = sprintf("select distinct(tag) from datasets where ");
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
    $sqlcmd = sprintf("%s and tag is not null order by tag",$sqlcmd);
    my $gettags = $dbh->prepare($sqlcmd);
    $gettags->execute();
    while (my @res = $gettags->fetchrow_array)
    {
	print "$res[0]\n";
    }
    exit(0);
}

sub printruns
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
    $sqlcmd = sprintf("%s and tag = '$tag' and dsttype='$ARGV[0]'order by runnumber",$sqlcmd);
    my $getruns =  $dbh->prepare($sqlcmd);
    $getruns->execute();
    if ($getruns->rows == 0)
    {
	print "no run found for tag $tag and dst type $ARGV[0]\n";
    }
    else
    {
	while (my @res = $getruns->fetchrow_array())
	{
	    print "$res[0]\n";
	}
    }
    exit(0);
}
