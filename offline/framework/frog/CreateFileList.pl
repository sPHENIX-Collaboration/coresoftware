#!/usr/bin/perl

#located in offline/framework/frog/

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(shuffle);

sub commonfiletypes;
sub fill_nocombine_files;
sub print_single_types;

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
#    "1" => "hijing (0-12fm) pileup 0-12fm DELETED",
#    "2" => "hijing (0-4.88fm) pileup 0-12fm DELETED",
    "3" => "pythia8 pp MB",
    "4" => "hijing (0-20fm) pileup 0-20fm",
#    "5" => "hijing (0-12fm) pileup 0-20fm DELETED",
    "6" => "hijing (0-4.88fm) pileup 0-20fm",
    "7" => "HF pythia8 Charm",
    "8" => "HF pythia8 Bottom",
    "9" => "HF pythia8 Charm D0",
    "10" => "HF pythia8 Bottom D0",
    "11" => "JS pythia8 Jet ptmin = 30GeV",
    "12" => "JS pythia8 Jet ptmin = 10GeV",
    "13" => "JS pythia8 Photon Jet",
    "14" => "Single Particles",
    "15" => "Special Productions",
    "16" => "HF pythia8 D0 Jets"
    );

my %pileupdesc = (
    "1" => "50kHz for Au+Au, 3MHz for p+p (default)",
    "2" => "25kHz for Au+Au",
    "3" => "10kHz for Au+Au"
    );

my $nEvents;
my $start_segment;
my $last_segment;
my $randomize;
my $prodtype;
my $runnumber = 40;
my $verbose;
my $nopileup;
my $embed;
my $pileup = 1;
my $particle;
my $pmin;
my $pmax;
my $production;

GetOptions('embed' => \$embed, 'l:i' => \$last_segment, 'n:i' => \$nEvents, "nopileup" => \$nopileup, "particle:s" => \$particle, 'pileup:i' => \$pileup, "pmin:i" => \$pmin, "pmax:i"=>\$pmax, "production:s"=>\$production, 'rand' => \$randomize, 'run:i' => \$runnumber, 's:i' => \$start_segment, 'type:i' =>\$prodtype, "verbose" =>\$verbose);
my $filenamestring;
my %filetypes = ();
my %notlike = ();

my $pileupstring;
my $pp_pileupstring;

if (defined $embed && defined $nopileup)
{
    print "--embed and --nopileup flags do not work together, it does not make sense\n";
    exit(1);
}
if ($pileup == 1)
{
    $pileupstring = sprintf("50kHz");
    $pp_pileupstring = sprintf("3MHz");
}
elsif ($pileup == 2)
{
    $pileupstring = sprintf("25kHz");
}
elsif ($pileup == 3)
{
    $pileupstring = sprintf("10kHz");
}
else
{
    print "invalid pileup option $pileup\n";
    exit(1);
}

my $embedok = 0;

if (defined $prodtype)
{
    if ($prodtype == 1)
    {
	$filenamestring = sprintf("sHijing_0_12fm_%s_bkg_0_12fm",$pileupstring);
        die "This dataset has been deleted\n";
	&commonfiletypes();
    }
    elsif ($prodtype == 2)
    {
	$filenamestring = sprintf("sHijing_0_488fm_%s_bkg_0_12fm",$pileupstring);
        die "Dataset $prodtype has been deleted\n";
	&commonfiletypes();
    }
    elsif ($prodtype == 3)
    {
	$filenamestring = "pythia8_pp_mb";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 4)
    {
	if (defined $nopileup)
	{
	    $filenamestring = sprintf("sHijing_0_20fm");
	}
	else
	{
	    $filenamestring = sprintf("sHijing_0_20fm_%s_bkg_0_20fm",$pileupstring);
	}
        $notlike{$filenamestring} = "pythia8";
	&commonfiletypes();
    }
    elsif ($prodtype == 5)
    {
	$filenamestring = sprintf("sHijing_0_12fm_%s_bkg_0_20fm",$pileupstring);
        die "Dataset $prodtype has been deleted\n";
	&commonfiletypes();
    }
    elsif ($prodtype == 6)
    {
	$filenamestring = sprintf("sHijing_0_488fm_%s_bkg_0_20fm",$pileupstring);
	&commonfiletypes();
    }
    elsif ($prodtype == 7)
    {
	$filenamestring = "pythia8_Charm";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 8)
    {
	$filenamestring = "pythia8_Bottom";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 9)
    {
	$filenamestring = "pythia8_CharmD0";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 10)
    {
	$filenamestring = "pythia8_BottomD0";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 11)
    {
        $embedok = 1;
	$filenamestring = "pythia8_Jet30";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		$filenamestring = sprintf("%s_sHijing_0_20fm_%s_bkg_0_20fm",$filenamestring, $pileupstring);
	    }
	    else
	    {
		$filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	    }
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 12)
    {
        $embedok = 1;
	$filenamestring = "pythia8_Jet10";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		$filenamestring = sprintf("%s_sHijing_0_20fm_%s_bkg_0_20fm",$filenamestring, $pileupstring);
	    }
	    else
	    {
		$filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	    }
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 13)
    {
        $embedok = 1;
	$filenamestring = "pythia8_PhotonJet";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		$filenamestring = sprintf("%s_sHijing_0_20fm_%s_bkg_0_20fm",$filenamestring, $pileupstring);
	    }
	    else
	    {
		$filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	    }
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 14)
    {
        $embedok = 1;
        $nopileup = 1;
        my $bad = 0;
	$filenamestring = "single";
	if (!defined $particle)
	{
	    print "-particle: G4 particle name needs to be set for single particle sims\n";
	    $bad = 1;
	}
	if ($bad > 0)
	{
            print "\nExisting single particle sims, use:\n";
            print_single_types();
	    exit(1);
	}
	if (defined $pmin && defined $pmax)
	{
	    $filenamestring = sprintf("%s_%s_%d_%dMeV",$filenamestring, $particle, $pmin, $pmax);
	}
	else
	{
	    if (defined $embed)
	    {
		$filenamestring = sprintf("%s_%s_sHijing_0_20fm_50kHz_bkg_0_20fm",$filenamestring, $particle);
	    }
	    else
	    {
		$filenamestring = sprintf("%s_%s",$filenamestring, $particle, $pmin, $pmax);
	    }
	}
	&commonfiletypes();
    }
    elsif ($prodtype == 15)
    {
        $nopileup = 1;
        my $bad = 0;
	$filenamestring = "special";
	if (! defined $production)
	{
	    $bad = 1;
	}
	if ($bad > 0)
	{
            print "\nExisting special sims, use:\n";
            print_special_types();
	    exit(1);
	}
	$filenamestring = sprintf("%s_%s",$filenamestring, $production);
	&commonfiletypes();
    }
    elsif ($prodtype == 16)
    {
	$filenamestring = "pythia8_JetD0";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s_%s",$filenamestring,$pp_pileupstring);
	}
	&commonfiletypes();
    }
    else
    {
	print "no production type $prodtype\n";
	exit(1);
    }
    &fill_other_types();
}

if (defined $embed && ! $embedok)
{
    print "Embedding not implemented for type $prodtype\n";
    exit(1);
}

my $filenamestring_with_runnumber = sprintf("%s\-%010d-",$filenamestring,$runnumber);
if ($#ARGV < 0)
{
    if (! defined $prodtype)
    {
	print "usage: CreateFileLists.pl -type <production type> <filetypes>\n";
	print "parameters:\n";
	print "-embed : pp embedded into hijing (only for pp types)\n";
	print "-l     : last segment\n";
	print "-n     : <number of events>\n";
	print "-nopileup : without pileup\n";
	print "-rand  : randomize segments used\n";
	print "-run   : runnumber (default = $runnumber)\n";
	print "-s     : starting segment (remember first segment is 0)\n";
	print "\n-type  : production type\n";
	foreach my $pd (sort { $a <=> $b } keys %proddesc)
	{
	    print "    $pd : $proddesc{$pd}\n";
	}
	print "\n-pileup : pileup rate selection (default = $pileup)\n";
	foreach my $pd (sort { $a <=> $b } keys %pileupdesc)
	{
	    print "    $pd : $pileupdesc{$pd}\n";
	}
        print "\n Single particle mandatory options:\n";
        print "-particle : G4 particle name\n";
        print "-pmin : minimum momentum (in MeV/c)\n";
        print "-pmax : maximum momentum (in MeV/c)\n";

        print "\n Special production mandatory options:\n";
        print "-production : production name\n";

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
print "This Can Take a While (10 minutes depending on the amount of events and the number of file types you want)\n";
my $conds = sprintf("dsttype = ? and filename like \'\%%%s\%\'",$filenamestring_with_runnumber);

if (exists $notlike{$filenamestring})
{
    $conds = sprintf("%s and filename not like  \'\%%%s\%\'",$conds,$notlike{$filenamestring});
}
if (defined $start_segment)
{
    $conds = sprintf("%s and segment >= %d",$conds,$start_segment);
}
if (defined $last_segment)
{
    if (defined $start_segment)
    {
	if ($last_segment < $start_segment)
	{
	    print "last segment: (-l $last_segment) smaller than start segment: (-s $start_segment), I will try but you will not get anything\n";
	}
    }
    $conds = sprintf("%s and segment <= %d",$conds,$last_segment);
}
my $getfilesql = sprintf("select filename,segment,events from datasets where %s order by segment",$conds);
#my $getfilesql = sprintf("select filename,segment,events from datasets where %s ",$conds);

my %getfiles = ();
foreach  my $tp (keys %req_types)
{
    if ($tp eq "G4Hits")
    {
	if (defined $embed)
	{
	    print "Selecting G4Hits with -embed is not supported (and does not make sense)\n";
	    exit(1);
	}
	my $newfilenamestring;
	if (defined $nopileup) # for no pileup we have the string already (G4Hits are nopileup)
	{
	    my @sp1 = split(/-/,$filenamestring_with_runnumber);
	    $newfilenamestring = $filenamestring_with_runnumber;
	}
	else
	{
	    my @sp1 = split(/_/,$filenamestring_with_runnumber);
	    if ($#sp1 == 3 || $#sp1 == 6 )
	    {
		$newfilenamestring = sprintf("%s_%s_%s\-%010d-",$sp1[0],$sp1[1],$sp1[2],$runnumber);
	    }
	    elsif ($#sp1 == 2)
	    {
		$newfilenamestring = sprintf("%s_%s\-%010d-",$sp1[0],$sp1[1],$runnumber);
	    }
	    else
	    {
		print "splitting $filenamestring_with_runnumber gave bad number of _: $#sp1\n";
		die;
	    }
	}
	my $newgetfilesql = $getfilesql;
	$newgetfilesql =~ s/$filenamestring_with_runnumber/$newfilenamestring/;
	$getfiles{"G4Hits"} = $dbh->prepare($newgetfilesql);
	if (defined $verbose)
	{
	    print "sql: $newgetfilesql\n";
	}
    }
    else
    {
	$getfiles{$tp} = $dbh->prepare($getfilesql);
	if (defined $verbose)
	{
	    print "sql: $getfilesql\n";
	}
    }
}
#die;
# here we fill the big hash with all segments/files for all requested filetypes
if (defined $verbose)
{
    print "fetching files from DB done, hashing all of them\n";
}
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
if (defined $verbose)
{
    print "hashing done, finding hash with lowest number of entries\n";
}
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
if (defined $verbose)
{
    print "matching hashes\n";
}

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

sub commonfiletypes
{
# pass1
    $filetypes{"G4Hits"} = "G4 Hits";
# pass2
    $filetypes{"DST_BBC_G4HIT"} = "Pileup BBC/MBD G4Hits";
    $filetypes{"DST_CALO_G4HIT"} = "Pileup Calorimeter G4Hits";
    $filetypes{"DST_TRKR_G4HIT"} = "Pileup Tracking Detector G4 Hits";
    $filetypes{"DST_TRUTH_G4HIT"} = "temporary Pileup Truth info, use DST_TRUTH";
    $filetypes{"DST_VERTEX"} = "Pileup Simulated Smeared Vertex";
# pass3 calo
    $filetypes{"DST_CALO_CLUSTER"} = "Reconstructed Calorimeter Towers and Clusters";
#pass3 trk
    $filetypes{"DST_TRKR_HIT"} = "TPC and Silicon Hits";
    $filetypes{"DST_TRUTH"} = "Truth Info (updated with Clusters)";
#pass4 tracks
    $filetypes{"DST_TRKR_CLUSTER"} = "pass0 output: tpc clusters";
    $filetypes{"DST_TRACKSEEDS"} = "passA output: track seeds";
    $filetypes{"DST_TRACKS"} = "passC output: Reconstructed Tracks";
#analysis pass
    $filetypes{"DST_TRUTH_JET"} = "Truth Jets";
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

sub print_single_types
{
    my $sqlstring = sprintf("select filename from datasets where runnumber = %d and filename like '%%_single_%%' and segment=0",$runnumber);
    my $getallfiles = $dbh->prepare($sqlstring);
    $getallfiles->execute();
    my $runsplit = sprintf("MeV-%010d",$runnumber);
    my $runsplit_embed = sprintf("_sHijing_0_20fm");
    my $runsplit_runnumber = sprintf("-%010d",$runnumber);
    my %types = ();
    my %dsts = ();
    while (my @res = $getallfiles->fetchrow_array())
    {
	my @sp1 = split(/_single_/,$res[0]);
        my @sp2;
        my $typeflag = "";
	if ($sp1[1] =~ /MeV/)
	{
	    @sp2 = split(/$runsplit/,$sp1[1]);
	}
	elsif ($sp1[1] =~ /$runsplit_embed/)
	{
	    @sp2 = split(/$runsplit_embed/,$sp1[1]);
            $typeflag = "-embed ";
	}
	else
	{
	    @sp2 = split(/$runsplit_runnumber/,$sp1[1]);
	}

	$types{$sp2[0]} = $typeflag;
        $dsts{$sp1[0]} = 1;
    }
    $getallfiles->finish();
    foreach my $name (sort keys %types)
    {
	if ($name =~ /(\S+)\_(\d+)\_(\d+).*/ )
	{
	    print "CreateFileList.pl -type 14 -run $runnumber -particle $1 -pmin $2 -pmax $3\n";
	}
	else
	{
	    print "CreateFileList.pl -type 14 $types{$name}-run $runnumber -particle $name\n";
	}
    }
    print "\nDST types:\n";
    foreach my $name (sort keys %dsts)
    {
	    print "$name\n";
    }
}

sub print_special_types
{
    my $sqlstring = sprintf("select filename from datasets where runnumber = %d and filename like '%%_special_%%'",$runnumber);
    my $getallfiles = $dbh->prepare($sqlstring);
    $getallfiles->execute();
    my $runsplit = sprintf("-%010d",$runnumber);
    my %types = ();
    my %dsts = ();
    while (my @res = $getallfiles->fetchrow_array())
    {
	my @sp1 = split(/_special_/,$res[0]);
	my @sp2 = split(/$runsplit/,$sp1[1]);
	$types{$sp2[0]} = 1;
        $dsts{$sp1[0]} = 1;
    }
    $getallfiles->finish();
    foreach my $name (sort keys %types)
    {
	    print "CreateFileList.pl -type 15 -production $name\n";
    }
    print "\nDST types:\n";
    foreach my $name (sort keys %dsts)
    {
	    print "$name\n";
    }
}
