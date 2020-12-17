#!/usr/bin/perl

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;

my $dbh = DBI->connect("dbi:ODBC:FileCatalog","argouser") || die $DBI::error;
$dbh->{LongReadLen}=2000; # full file paths need to fit in here

my $getdsttypes = $dbh->prepare("select distinct(dsttype) from datasets");
$getdsttypes->execute();

my %dsttype = ();
while(my @res = $getdsttypes->fetchrow_array())
{
    my $listfile = sprintf("%s.list",lc $res[0]);
    $dsttype{$res[0]} = $listfile;
    if (-f $listfile)
    {
	unlink $listfile;
    }
}

if ($#ARGV < 0)
{
    print "usage: CreateFileLists.pl -fm <fermi range> <filetypes>\n";
    print "parameters:\n";
    print "-fm : fermi range, 0-12 = 1, 0-4.88 = 2\n";
    print "-n : <number of events>\n";
    print "-s : <starting segment>\n";
#    print "-r : pick random segments until nEvents is reached\n";
    print "\navailable types:\n";
    foreach my $tp (sort keys %dsttype)
    {
	print "$tp\n";
    }
    exit(0);
}
my $nEvents;
my $start_segment;
my $randomize;
my $fmrange;
GetOptions('fm:i' =>\$fmrange, 'n=i' => \$nEvents, 'r' => \$randomize, 's=i' => \$start_segment);

if (! defined $fmrange || $fmrange < 1 || $fmrange > 1)
{
    if (defined $fmrange)
    {
	print "-fm $fmrange invalid, need to give valid fermi range\n";
    }
    else
    {
	print "need to give fermi range\n";
    }
    print "-fm 1 : 0-12fm\n";
    print "-fm 2 : 0-4.88fm\n";
    exit(0);
}
my $fmstring;
if ($fmrange == 1)
{
    $fmstring = "0_12fm";
}
elsif ($fmrange == 2)
{
    $fmstring = "0_488fm";
}
else
{
    print "no string for fermi range $fmrange\n";
    exit(1);
}
#check dst types
my %req_types = ();
my %allfilehash = ();
my %allevthash = ();

while($#ARGV >= 0)
{
    if (! exists $dsttype{$ARGV[0]})
    {
	print "dst type $ARGV[0] does not exist\n";
	print "\navailable types:\n";
	foreach my $tp (sort keys %dsttype)
	{
	    print "$tp\n";
	}
	exit(1);
    }
    $req_types{$ARGV[0]} = 1;
    $allfilehash{$ARGV[0]} = ();
    $allevthash{$ARGV[0]} = ();
    shift (@ARGV);
}
my $conds = sprintf("dsttype = ? and filename like '%$fmstring%'");
if (defined $start_segment)
{
    $conds = sprintf("%s and segment >= %d",$conds,$start_segment);
}
my $getfilesql = sprintf("select filename,segment,events from datasets where %s order by segment",$conds);

#my $getfiles = $dbh->prepare("select filename,segment,events from datasets where dsttype = ? and filename like '%$fmstring%' order by segment");
my $getfiles = $dbh->prepare($getfilesql);

foreach my $tp (sort keys %req_types)
{
    my %dsthash = ();
    my %evthash = ();
    $getfiles->execute($tp);
    if ($getfiles->rows == 0)
    {
	print "no files for type $tp\n";
	exit(0);
    }
    while (my @res = $getfiles->fetchrow_array())
    {
	my $hashkey = sprintf("%05d",$res[1]);
	$dsthash{$hashkey} = $res[0];
	$evthash{$res[0]} = $res[2];
    }
    $allfilehash{$tp} = \%dsthash;
    $allevthash{$tp} = \%evthash;
}

my $entries = 100000;
my $lowtype;
foreach my $tp (sort keys %allfilehash)
{
    if ($entries > keys %{$allfilehash{$tp}})
    {
	$entries = keys %{$allfilehash{$tp}};
	$lowtype = $tp;
    }
}
#print "lowest entries: $entries, type: $lowtype\n";
#print Dumper(%allevthash);
foreach my $seg (sort keys %{$allfilehash{$lowtype}})
{
    my $isgood = 1;
    foreach my $tp (sort keys %allfilehash)
    {
	if ($tp eq $lowtype)
	{
	    next;
	}
	if (! exists  $allfilehash{$tp}{$seg})
	{
	    $isgood = 0;
	    last;
	}
    }

    if ($isgood == 1)
    {
#	print "segment $seg is good\n";
	foreach my $tp (sort keys %allfilehash)
	{
#	    print "using $allfilehash{$tp}{$seg}\n";
	    my $printcmd = sprintf("echo %s >> %s",$allfilehash{$tp}{$seg},$dsttype{$tp});
	    system($printcmd);
	}

	if (defined $nEvents)
	{
#	print "events: from $allfilehash{$lowtype}{$seg}: $allevthash{$lowtype}{$allfilehash{$lowtype}{$seg}}\n";
	    $nEvents = $nEvents - $allevthash{$lowtype}{$allfilehash{$lowtype}{$seg}};
	    if ($nEvents <= 0)
	    {
		last;
	    }
	}
    }
    else
    {
	print "segment $seg is bad\n";
    }

}


$getdsttypes->finish();
$getfiles->finish();
$dbh->disconnect;
