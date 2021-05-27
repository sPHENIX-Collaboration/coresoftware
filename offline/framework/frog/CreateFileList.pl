#!/usr/bin/perl

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(shuffle);

sub hijingfiletypes;
sub charmfiletypes;
sub bottomfiletypes;
sub pythiambfiletypes;
sub fill_nocombine_files;

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

my %proddesc = (
    "1" => "hijing (0-12fm) pileup 0-12fm",
    "2" => "hijing (0-4.88fm) pileup 0-12fm",
    "3" => "pythia8 pp MB",
    "4" => "hijing (0-20fm) pileup 0-20fm",
    "5" => "hijing (0-12fm) pileup 0-20fm",
    "6" => "hijing (0-4.88fm) pileup 0-20fm",
    "7" => "HF pythia8 Charm",
    "8" => "HF pythia8 Bottom"
    );


my $nEvents;
my $start_segment;
my $randomize;
my $prodtype;
GetOptions('type:i' =>\$prodtype, 'n:i' => \$nEvents, 'r' => \$randomize, 's:i' => \$start_segment);

my $filenamestring;
my %filetypes = ();
if (defined $prodtype)
{
    if ($prodtype == 1)
    {
	$filenamestring = "sHijing_0_12fm_50kHz_bkg_0_12fm-";
	&hijingfiletypes();
    }
    elsif ($prodtype == 2)
    {
	$filenamestring = "sHijing_0_488fm_50kHz_bkg_0_12fm-";
	&hijingfiletypes();
    }
    elsif ($prodtype == 3)
    {
	$filenamestring = "pythia8_mb";
	&pythiambfiletypes();
    }
    elsif ($prodtype == 4)
    {
	$filenamestring = "sHijing_0_20fm_50kHz_bkg_0_20fm-";
	&hijingfiletypes();
    }
    elsif ($prodtype == 5)
    {
	$filenamestring = "sHijing_0_12fm_50kHz_bkg_0_20fm-";
	&hijingfiletypes();
    }
    elsif ($prodtype == 6)
    {
	$filenamestring = "sHijing_0_488fm_50kHz_bkg_0_20fm-";
	&hijingfiletypes();
    }
    elsif ($prodtype == 7)
    {
	$filenamestring = "DST_HF_CHARM";
	&charmfiletypes();
    }
    elsif ($prodtype == 8)
    {
	$filenamestring = "DST_HF_BOTTOM";
	&bottomfiletypes();
    }
    else
    {
	print "no file substring for production type $prodtype\n";
	exit(1);
    }
    &fill_other_types();
}

if ($#ARGV < 0)
{
    if (! defined $prodtype)
    {
	print "usage: CreateFileLists.pl -type <production type> <filetypes>\n";
	print "parameters:\n";
	print "-n  : <number of events>\n";
	print "-r  : randomize segments used\n";
	print "-s  : <starting segment>\n";
	print "-type : production type\n";
	foreach my $pd (sort keys %proddesc)
	{
	    print "    $pd : $proddesc{$pd}\n";
	}
	print "\navailable file types (choose at least one, --> means: written to):\n";
	foreach my $tp (sort keys %dsttype)
	{
	    print "$tp  --> $dsttype{$tp}\n";
	}
    }
    else
    {
	print "\navailable file types for -type $prodtype: $proddesc{$prodtype}:\n";
	foreach my $tp (sort keys %filetypes)
	{
	    print "\t$tp : $filetypes{$tp}\n";
	}
    }

    exit(0);
}

if (defined $randomize && ! defined $nEvents)
{
    print "randomizing segments only if number of events is selected\n";
    exit(0);
}

if (! defined $prodtype )
{
    print "need to give production type\n";
    print "-type : production type\n";
    foreach my $pd (sort keys %proddesc)
    {
	print "    $pd : $proddesc{$pd}\n";
    }
    exit(0);
}
if (! exists $proddesc{$prodtype})
{
    print  "invalid production type $prodtype, valid values\n";
    print "-type : production type\n";
    foreach my $pd (sort keys %proddesc)
    {
	print "    $pd : $proddesc{$pd}\n";
    }
    exit(0);
}

#check dst types
my %req_types = ();

# hash of {dsttype} containing hash of segment,filename
my %allfilehash = ();

my %allevthash = ();

my %nocombine = ();
&fill_nocombine_files;

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
    if (exists $req_types{$ARGV[0]})
    {
	print "please no duplicate file types ($ARGV[0])\n";
	exit(1);
    }

    if (exists $nocombine{$ARGV[0]})
    {
	if ($#ARGV >= 1 || keys %req_types > 0)
	{
	    print "File type $ARGV[0] cannot be combined with other files\n";
	    exit(1);
	}
    }

    $req_types{$ARGV[0]} = 1;
    $allfilehash{$ARGV[0]} = ();
    $allevthash{$ARGV[0]} = ();
    shift (@ARGV);

}
print "This Can Take a While (a minute give or take)\n";
my $conds = sprintf("dsttype = ? and filename like \'\%%%s\%\'",$filenamestring);
if (defined $start_segment)
{
    $conds = sprintf("%s and segment >= %d",$conds,$start_segment);
}
my $getfilesql = sprintf("select filename,segment,events from datasets where %s order by segment",$conds);

#print "sql: $getfilesql\n";

my %getfiles = ();
foreach  my $tp (keys %req_types)
{
    if ($tp eq "G4Hits")
    {
	my @sp1 = split(/_/,$filenamestring);
	my $newfilenamestring = sprintf("%s_%s_%s",$sp1[0],$sp1[1],$sp1[2]);
	my $newgetfilesql = $getfilesql;
	$newgetfilesql =~ s/$filenamestring/$newfilenamestring/;
	$getfiles{"G4Hits"} = $dbh->prepare($newgetfilesql);
    }
    else
    {
	$getfiles{$tp} = $dbh->prepare($getfilesql);
    }
}
#die;
# here we fill the big hash with all segments/files for all requested filetypes
foreach my $tp (sort keys %req_types)
{
    my %dsthash = ();
    my %evthash = ();
    $getfiles{$tp}->execute($tp);
    if ($getfiles{$tp}->rows == 0)
    {
	print "no files for type $tp\n";
	exit(0);
    }
    while (my @res = $getfiles{$tp}->fetchrow_array())
    {
	my $hashkey = sprintf("%05d",$res[1]);
	$dsthash{$hashkey} = $res[0];
	$evthash{$res[0]} = $res[2];
    }
    $allfilehash{$tp} = \%dsthash;
    $allevthash{$tp} = \%evthash;
}

my $entries = 100000; # given that we have 1000 files max, this value is always higher
my $lowtype;
# here we find the dst type with the smallest number of entries (segments)
# so we do not loop too much when finding matches for the other types
foreach my $tp (sort keys %allfilehash)
{
    if ($entries > keys %{$allfilehash{$tp}})
    {
	$entries = keys %{$allfilehash{$tp}};
	$lowtype = $tp;
    }
}
# here $lowtype is the dst type with the smallest number of segments
#print "lowest entries: $entries, type: $lowtype\n";
#print Dumper(%allevthash);

my @segarray = ();
foreach my $seg (sort keys %{$allfilehash{$lowtype}})
{
    foreach my $tp (sort keys %allfilehash)
    {
	if ($tp eq $lowtype)
	{
	    next;
	}
	if (! exists  $allfilehash{$tp}{$seg})
	{
	    last;
	}
    }
    push(@segarray,$seg);
}
# in segarray we have the common segments of all files now

# remove segments from array when number of events is reached
# randomize segments when -r is set
if (defined $nEvents)
{
    if (defined $randomize)
    {
	@segarray = shuffle(@segarray);
    }
    my @tmparray = ();
    foreach my $seg (@segarray)
    {
	push(@tmparray,$seg);
	$nEvents -= $allevthash{$lowtype}{$allfilehash{$lowtype}{$seg}};
	if ($nEvents <= 0)
	{
	    last;
	}
    }
    @segarray = @tmparray;
}
# sort list of segments and write to output file
my $nSelectedEvents = 0;
foreach my $seg (sort @segarray)
{
    $nSelectedEvents += $allevthash{$lowtype}{$allfilehash{$lowtype}{$seg}};
#	print "segment $seg is good\n";
    foreach my $tp (sort keys %allfilehash)
    {
#	    print "using $allfilehash{$tp}{$seg}\n";
	my $printcmd = sprintf("echo %s >> %s",$allfilehash{$tp}{$seg},$dsttype{$tp});
	system($printcmd);
    }

}
print "wrote the following list files containing >= $nSelectedEvents events:\n";
foreach my $tp (sort keys %allfilehash)
{
    print "$dsttype{$tp}\n";
}


$getdsttypes->finish();
foreach my $tp (keys %getfiles)
{
    $getfiles{$tp}->finish();
}
$dbh->disconnect;

sub hijingfiletypes
{
# pass1
    $filetypes{"G4Hits"} = "G4 Hits";
# pass2
    $filetypes{"DST_BBC_G4HIT"} = "Pileup BBC/MBD G4Hits";
    $filetypes{"DST_CALO_G4HIT"} = "Pileup Calorimeter G4Hits";
    $filetypes{"DST_TRKR_G4HIT"} = "Pileup Tracking Detector G4 Hits";
    $filetypes{"DST_TRUTH_G4HIT"} = "Pileup Truth info";
    $filetypes{"DST_VERTEX"} = "Pileup Simulated Smeared Vertex";
# pass3 calo
    $filetypes{"DST_CALO_CLUSTER"} = "Reconstructed Calorimeter Towers and Clusters";
#pass3 trk
    $filetypes{"DST_TRKR_CLUSTER"} = "TPC/Silicon Clusters";
#pass4 tracks
    $filetypes{"DST_TRACKS"} = "Reconstructed Tracks";
}

sub charmfiletypes
{
    $filetypes{"DST_HF_CHARM"} = "Charm DST";
    $filetypes{"QA_DST_HF_CHARM"} = "Charm QA";
    $filetypes{"JET_EVAL_DST_HF_CHARM"} = "Charm Jet Eval";
}

sub bottomfiletypes
{
    $filetypes{"DST_HF_BOTTOM"} = "Bottom DST";
    $filetypes{"QA_DST_HF_BOTTOM"} = "Bottom QA";
    $filetypes{"JET_EVAL_DST_HF_BOTTOM"} = "Bottom Jet Eval";
}

sub pythiambfiletypes
{
    $filetypes{"G4Hits"} = "G4 Hits";
# pass2
    $filetypes{"DST_BBC_G4HIT"} = "Pileup BBC/MBD G4Hits";
    $filetypes{"DST_CALO_G4HIT"} = "Pileup Calorimeter G4Hits";
    $filetypes{"DST_TRKR_G4HIT"} = "Pileup Tracking Detector G4 Hits";
    $filetypes{"DST_TRUTH_G4HIT"} = "Pileup Truth info";
    $filetypes{"DST_VERTEX"} = "Pileup Simulated Smeared Vertex";
}

# here are filetypes which are ntuples or cannot be combined
# with other files for any other reason
sub fill_nocombine_files
{
    $nocombine{"JET_EVAL_DST_HF_CHARM"} = 1;
    $nocombine{"JET_EVAL_DST_HF_BOTTOM"} = 1;
    $nocombine{"QA_DST_HF_CHARM"} = 1;
    $nocombine{"QA_DST_HF_BOTTOM"} = 1;
}

sub fill_other_types
{
    my $sqlstring = sprintf("select distinct(dsttype) from datasets where filename like '%%%s%%'",$filenamestring);
    my $getalltypes = $dbh->prepare($sqlstring);
    $getalltypes->execute();
    while (my @res = $getalltypes->fetchrow_array())
    {
	if (! exists $filetypes{$res[0]})
	{
	    $filetypes{$res[0]} = "No Description";
	}
    }
    $getalltypes->finish();
}
