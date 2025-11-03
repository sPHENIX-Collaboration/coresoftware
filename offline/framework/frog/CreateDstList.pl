#!/usr/bin/env perl

#located in offline/framework/frog/

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;

sub creatfilelists;
sub printdatasets;
sub printdsttypes;
sub printtags;
sub printruns;

my $tag;
my $verbose;
my $runnumber;
my $runlist;
my $printdatasets;
my $printdsttypes;
my $printruns;
my $printtags;
my $hpss;
my $noprint;
my $dataset;

GetOptions('dataset:s' => \$dataset, 'hpss' =>\$hpss, 'list:s' => \$runlist, 'printdatasets' => \$printdatasets, 'printdsttypes' => \$printdsttypes, 'printruns' => \$printruns, 'printtags' => \$printtags, 'run:i' => \$runnumber, 'tag:s' => \$tag, "verbose" =>\$verbose);

my %ignore_datasets = (
    "mdc2" => 1,
    );

my $dbh = DBI->connect("dbi:ODBC:FileCatalog_read") || die $DBI::error;
$dbh->{LongReadLen}=2000; # full file paths need to fit in here

if (defined $printdatasets)
{
    printdatasets();
}

if (defined $printdsttypes)
{
    printdsttypes();
}

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
	print "--hpss <print files in hpss>\n";
	print "--list <file with list of runs>\n";
	print "--printdatasets: print existing datasets\n";
	print "--printdsttypes: print existing dst types\n";
        print "             needs --dataset <dataset> flag\n";
	print "--printruns: print existing runs (for piping into a runlist)\n";
        print "             needs --dataset <dataset> flag\n";
        print "             needs --tag <tag> flag\n";
	print "--printtags: print existing tags\n";
	print "--run <run number>\n";
	print "--tag <tag (build_cdb_version)> (mandatory)\n";
	print "--version <version> (only if we have multiple versions for same build/cdb tag)\n";
	print "--verbose print stuff\n";
	print "\nfor use with --printruns or --printtags\n";
        print "--dataset <dataset (run2pp, run2auau, run3auau, run3cosmics,...)\n"; 
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
	else
	{
	    print "wrote $filename\n";
	}
    }
    #x	print Dumper(\%files);
}

sub printtags
{
    my $sqlcmd = sprintf("select datasets.tag,datasets.dataset from datasets group by tag,dataset");
    my $icnt = 0;
    $sqlcmd = sprintf("%s",$sqlcmd);
    if (defined $verbose)
    {
	print "sql cmd: $sqlcmd\n";
    }
    my $gettags = $dbh->prepare($sqlcmd);
    $gettags->execute();
    my %tags = ();
    while (my @res = $gettags->fetchrow_array)
    {
	if (defined $dataset)
	{
	    if ($res[1] eq $dataset)
	    {
		$tags{$res[0]} = $res[1];
	    }
	}
	else
	{
	    my $addthis = 1;
	    foreach my $ignore (keys %ignore_datasets)
	    {
		if ($res[1] eq $ignore)
		{
		    $addthis = 0;
		    last;
		}
	    }
	    if (! defined $res[0])
	    {
		next;
	    }
	    if ($addthis == 1)
	    {
		$tags{$res[0]} = $res[1];
	    }
	}
    }
    foreach my $tg (sort keys %tags)
    {
	print "$tg\n";
    }
    exit(0);
}

sub printruns
{
    my %dsttype = ();
    while($#ARGV >= 0)
    {
	$dsttype{$ARGV[0]} = 1;
	shift (@ARGV);
    }
    my $isgood = 1;
    if (!defined $dataset)
    {
	print "printruns needs --dataset <dataset>\n";
	$isgood = 0;
    }
    if (! defined $tag)
    {
	print "printruns needs --tag <tag>\n";
	$isgood = 0;
    }
    if ($isgood == 0)
    {
	exit 1;
    }
    my $sqlcmd = sprintf("select datasets.runnumber,datasets.dsttype from datasets where dataset = '%s' and tag = '%s' group by runnumber,dsttype",$dataset,$tag);
    my $getruns =  $dbh->prepare($sqlcmd);
    $getruns->execute();
    my %selruns = ();
    if ($getruns->rows == 0)
    {
	print "no run found for tag $tag and dataset $dataset\n";
    }
    else
    {
	my $icnt = 0;
	while (my @res = $getruns->fetchrow_array())
	{
	    if (keys %dsttype > 0 && !exists $dsttype{$res[1]})
	    {
		next;
	    }
	    $selruns{$res[0]} = 1;
	    $icnt++;
	}
	if ($icnt == 0)
	{
	    print "no run for selected dst types\n";
	    exit 1;
	}
	foreach my $rn (sort keys %selruns)
	{
	    print "$rn\n";
	}
    }
    exit(0);
}

sub printdatasets
{
    my $sqlcmd = sprintf("select datasets.dataset from datasets group by dataset");
    my $getdatasets = $dbh->prepare($sqlcmd);
    $getdatasets->execute();
    my %ds = ();
    while (my @res = $getdatasets->fetchrow_array())
    {
	my $addthis = 1;
	foreach my $ignore (keys %ignore_datasets)
	{
	    if ($res[0] eq $ignore)
	    {
		$addthis = 0;
		last;
	    }
	}
	if ($addthis == 1)
	{
	    $ds{$res[0]} = 1;
	}
    }
    foreach my $dd (sort keys %ds)
    {
	print "$dd\n";
    }
    exit(0);
}

sub printdsttypes
{
    if (!defined $dataset)
    {
	print "printruns needs dataset\n";
    }
    my $sqlcmd = sprintf("select datasets.dsttype,datasets.tag from datasets where dataset = '%s' group by dsttype,tag",$dataset);
    my $getdsttypes = $dbh->prepare($sqlcmd);
    $getdsttypes->execute();
    my %dsts = ();
    while (my @res = $getdsttypes->fetchrow_array())
    {
	if (defined $tag)
	{
	    if ($res[1] ne $tag)
	    {
		next;
	    }
	}
	$dsts{$res[0]} = 1;
    }
    foreach my $dd (sort keys %dsts)
    {
	print "$dd\n";
    }
    exit(0);
}
