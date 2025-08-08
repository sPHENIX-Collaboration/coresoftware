#!/usr/bin/env perl

#located in offline/framework/frog/

use DBI;
use strict;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(shuffle);

sub commonfiletypes;
sub fill_nocombine_files;
sub print_single_types;
sub print_runs;

my $dbh = DBI->connect("dbi:ODBC:FileCatalog_read") || die $DBI::error;
$dbh->{LongReadLen}=2000; # full file paths need to fit in here

my $getdsttypes = $dbh->prepare("select distinct(dsttype) from datasets where dsttype not like '%\_pi\_%' ESCAPE '\' and dsttype <> 'beam' and dsttype <> 'cosmic' and dataset = 'mdc2'");
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
my %exclude_these = (
    "DST_JOBA" => "Test PanDa",
    "DST_MDC2_GLOBAL" => "Test PanDa",
    "DST_PASS1_CLUSTERS" => "Test PanDa",
    "DST_RECO_CLUSTER" => "Test PanDa"
    );

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
    "16" => "HF pythia8 D0 Jets",
    "17" => "HF pythia8 D0 pi-k Jets ptmin = 5GeV ",
    "18" => "HF pythia8 D0 pi-k Jets ptmin = 12GeV",
    "19" => "JS pythia8 Jet ptmin = 40GeV",
    "20" => "hijing pAu (0-10fm) pileup 0-10fm",
    "21" => "JS pythia8 Jet ptmin = 20GeV",
    "22" => "cosmic field on",
    "23" => "cosmic field off",
    "24" => "AMPT",
    "25" => "EPOS",
    "26" => "JS pythia8 Detroit",
    "27" => "JS pythia8 Photonjet ptmin = 5GeV",
    "28" => "JS pythia8 Photonjet ptmin = 10GeV",
    "29" => "JS pythia8 Photonjet ptmin = 20GeV",
    "30" => "Herwig MB",
    "31" => "Herwig Jet ptmin = 10 GeV",
    "32" => "Herwig Jet ptmin = 30 GeV",
    "33" => "JS pythia8 Jet ptmin = 15GeV",
    "34" => "JS pythia8 Jet ptmin = 50GeV"
    );

my %pileupdesc = (
    "1" => "50kHz for Au+Au, 3MHz for p+p (default)",
    "2" => "25kHz for Au+Au",
    "3" => "10kHz for Au+Au",
    "4" => "1MHz for pp 100us streaming",
    "5" => "2MHz for pp 20us streaming",
    ">5" => "pileup rate in kHz"
    );

my $nEvents;
my $start_segment;
my $last_segment;
my $randomize;
my $prodtype;
my $runnumber;
my $verbose;
my $nopileup;
my $nobkgpileup;
my $embed;
my $pileup = 1;
my $particle;
my $pmin;
my $pmax;
my $production;
my $momentum;
# that should teach me a lesson to not give a flag an optional string value
# just using embed:s leads to the next ARGV to be used as argument, even if it
# is the next option. Sadly getopt swallows the - so parsing this becomes
# quickly a nightmare, The only solution I see is to read the ARGV's - check for
# the argument in question with a string compare (=~/em/) and then have a look
# at the next argument (if it exists) and check if is is another option (=~/-/)
# and if so use "auau" as default for $embed and push the modified option
# into the new command line ARGV. It defaults to auau if -emb is set but
# neither pau nor auau is given
my @newargs = ();
my $iarg = 0;
foreach my $argument (@ARGV)
{
    if (substr($argument,1,2) eq "em")
    {
	my $firstchar = substr($ARGV[$iarg+1],0,1);
	if (! exists $ARGV[$iarg+1] || substr($ARGV[$iarg+1],0,1) eq "-")
	{
	    push(@newargs, $argument);
	    push(@newargs,"auau");
	}
	else
	{
	    push(@newargs, $argument);
	    if ($ARGV[$iarg+1] ne "pau" && $ARGV[$iarg+1] ne "auau" && $ARGV[$iarg+1] ne "central")
	    {
		push(@newargs,"auau");
	    }
	}
    }
    else
    {
	push(@newargs,$argument);
    }
    $iarg++;
}
@ARGV=@newargs;
GetOptions('embed:s' => \$embed, 'l:i' => \$last_segment, 'momentum:s' => \$momentum, 'n:i' => \$nEvents, "nobkgpileup" => \$nobkgpileup, "nopileup" => \$nopileup, "particle:s" => \$particle, 'pileup:i' => \$pileup, "pmin:i" => \$pmin, "pmax:i"=>\$pmax, "production:s"=>\$production, 'rand' => \$randomize, 'run:i' => \$runnumber, 's:i' => \$start_segment, 'type:i' =>\$prodtype, "verbose" =>\$verbose);
my $filenamestring;
my %filetypes = ();
my %notlike = ();

my $AuAu_pileupstring;
my $pp_pileupstring;
my $pAu_pileupstring;
my $pileupstring;

if (! defined $runnumber && $#newargs >= 0)
{
    print "\nyou need to give a runnumber with -run <runnumber>\n";
    print_runs();
    exit(1);
}
if (defined $embed && defined $nopileup)
{
    print "--embed and --nopileup flags do not work together, it does not make sense\n";
    exit(1);
}
if (!defined $embed && defined $nobkgpileup)
{
    print "--nobkgpileup flag only valid for embedding (use also --embed)\n";
    exit(1);
}
my $pAu_bkgpileup = sprintf("_bkg_0_20fm");
my $AuAu_bkgpileup = sprintf("_bkg_0_10fm");
if ($pileup == 1)
{
    $AuAu_pileupstring = sprintf("_50kHz%s",$AuAu_bkgpileup);
    $pp_pileupstring = sprintf("_3MHz");
    $pAu_pileupstring = sprintf("_500kHz%s",$pAu_bkgpileup);
}
elsif ($pileup == 2)
{
    $AuAu_pileupstring = sprintf("_25kHz%s",$AuAu_bkgpileup);
}
elsif ($pileup == 3)
{
    $AuAu_pileupstring = sprintf("_10kHz%s",$AuAu_bkgpileup);
}
elsif ($pileup == 4)
{
    $pp_pileupstring = sprintf("_1MHz");
}
elsif ($pileup == 5)
{
    $pp_pileupstring = sprintf("_2MHz");
}
else
{
    $pp_pileupstring = sprintf("_%dkHz",$pileup);
    $AuAu_pileupstring = sprintf("_%dkHz%s",$AuAu_bkgpileup);
}
if (defined $nobkgpileup)
{
    $pp_pileupstring = sprintf("");
    $AuAu_pileupstring = sprintf("");
}

my $embedok = 0;

if (defined $prodtype)
{
    if ($prodtype == 1)
    {
        die "This dataset has been deleted\n";
	&commonfiletypes();
    }
    elsif ($prodtype == 2)
    {
        die "Dataset $prodtype has been deleted\n";
	&commonfiletypes();
    }
    elsif ($prodtype == 3)
    {
	$filenamestring = "pythia8_pp_mb";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
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
	    $filenamestring = sprintf("sHijing_0_20fm%s",$AuAu_pileupstring);
	}
        $notlike{$filenamestring} = ["pythia8" ,"single", "special"];
        $pileupstring = $AuAu_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 5)
    {
	$filenamestring = sprintf("sHijing_0_12fm%s",$AuAu_pileupstring);
        die "Dataset $prodtype has been deleted\n";
	&commonfiletypes();
    }
    elsif ($prodtype == 6)
    {
	if (defined $nopileup)
	{
	    $filenamestring = sprintf("sHijing_0_488fm");
	}
	else
	{
	  $filenamestring = sprintf("sHijing_0_488fm%s",$AuAu_pileupstring);
	}
        $notlike{$filenamestring} = ["pythia8" ,"single", "special"];
        $pileupstring = $AuAu_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 7)
    {
	$filenamestring = "pythia8_Charm";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 8)
    {
	$filenamestring = "pythia8_Bottom";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 9)
    {
	$filenamestring = "pythia8_CharmD0";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 10)
    {
	$filenamestring = "pythia8_BottomD0";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
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
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
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
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
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
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
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
	    if (defined $momentum)
	    {
		$filenamestring = sprintf("%s_%s_%s_%d_%dMeV",$filenamestring, $particle, $momentum, $pmin, $pmax);
	    }
	    else
	    {
		$filenamestring = sprintf("%s_%s_%d_%dMeV",$filenamestring, $particle, $pmin, $pmax);
	    }

	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	}
	else
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
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
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 17)
    {
	$filenamestring = "pythia8_CharmD0piKJet5";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 18)
    {
	$filenamestring = "pythia8_CharmD0piKJet12";
	if (! defined $nopileup)
	{
	    $filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 19)
    {
        $embedok = 1;
	$filenamestring = "pythia8_Jet40";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 20)
    {
	if (defined $nopileup)
	{
	    $filenamestring = sprintf("sHijing_pAu_0_10fm");
	}
	else
	{
	    $filenamestring = sprintf("sHijing_pAu_0_10fm%s",$pAu_pileupstring);
	}
        $notlike{$filenamestring} = ["pythia8" ,"single", "special"];
        $pileupstring = $pAu_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 21)
    {
        $embedok = 1;
	$filenamestring = "pythia8_Jet20";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 22)
    {
	$filenamestring = "cosmic_magnet_on";
        $filenamestring = sprintf("%s",$filenamestring);
	$nopileup = 1; # it is no pileup only - no need to require it
	&commonfiletypes();
    }
    elsif ($prodtype == 23)
    {
	$filenamestring = "cosmic_magnet_off";
        $filenamestring = sprintf("%s",$filenamestring);
	$nopileup = 1; # it is no pileup only - no need to require it
	&commonfiletypes();
    }
    elsif ($prodtype == 24)
    {
	if (defined $nopileup)
	{
	    $filenamestring = sprintf("ampt_0_20fm");
	}
	else
	{
	    $filenamestring = sprintf("ampt_0_20fm%s",$AuAu_pileupstring);
	}
        $notlike{$filenamestring} = ["pythia8" ,"single", "special"];
        $pileupstring = $AuAu_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 25)
    {
	if (defined $nopileup)
	{
	    $filenamestring = sprintf("epos_0_153fm");
	}
	else
	{
	    $filenamestring = sprintf("epos_0_153fm_%s_bkg_0_153fm",$AuAu_pileupstring);
	}
        $notlike{$filenamestring} = ["pythia8" ,"single", "special"];
        $pileupstring = $pAu_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 26)
    {
        $embedok = 1;
	$filenamestring = "pythia8_Detroit";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 27)
    {
        $embedok = 1;
	$filenamestring = "pythia8_PhotonJet5";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 28)
    {
        $embedok = 1;
	$filenamestring = "pythia8_PhotonJet10";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 29)
    {
        $embedok = 1;
	$filenamestring = "pythia8_PhotonJet20";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 30)
    {
        $embedok = 1;
	$filenamestring = "Herwig_MB";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 31)
    {
        $embedok = 1;
	$filenamestring = "Herwig_Jet10";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 32)
    {
        $embedok = 1;
	$filenamestring = "Herwig_Jet30";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 33)
    {
        $embedok = 1;
	$filenamestring = "pythia8_Jet15";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
	&commonfiletypes();
    }
    elsif ($prodtype == 34)
    {
        $embedok = 1;
	$filenamestring = "pythia8_Jet50";
	if (! defined $nopileup)
	{
	    if (defined $embed)
	    {
		if ($embed eq "pau")
		{
		    $filenamestring = sprintf("%s_sHijing_pAu_0_10fm%s",$filenamestring, $pAu_pileupstring);
		}
		elsif ($embed eq "central")
		{
		    $filenamestring = sprintf("%s_sHijing_0_488fm%s",$filenamestring, $AuAu_pileupstring);
		}
		else
		{
		    $filenamestring = sprintf("%s_sHijing_0_20fm%s",$filenamestring, $AuAu_pileupstring);
		}
	    }
	    else
	    {
		$filenamestring = sprintf("%s%s",$filenamestring,$pp_pileupstring);
	    }
	}
        $pileupstring = $pp_pileupstring;
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
	print "-embed : pp embedded into MB AuAu hijing (only for pp types)\n";
	print "  -embed pau : embedded into pAu (only for pp types)\n";
	print "  -embed central : embedded into central AuAu\n";
	print "-l     : last segment\n";
	print "-n     : <number of events>\n";
	print "-nopileup : without pileup\n";
	print "-rand  : randomize segments used\n";
	print "-run   : runnumber (mandatory, no default anymore)\n";
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
	print "\n-nobkgpileup : no pileup of background event (use with -embed)\n";
        print "\n Single particle mandatory options:\n";
        print "-particle : G4 particle name\n";
        print "-mom : (optional) p or pt\n";
        print "-pmin : minimum momentum (in MeV/c)\n";
        print "-pmax : maximum momentum (in MeV/c)\n";

        print "\n Special production mandatory options:\n";
        print "-production : production name\n";

	print "\navailable file types (choose at least one, --> means: written to):\n";
	foreach my $tp (sort keys %dsttype)
	{
	    if (! exists $exclude_these{$tp})
	    {
		print "$tp  --> $dsttype{$tp}\n";
	    }
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
    my $ref = $notlike{$filenamestring};
    foreach my $item  (@$ref)
    {
	$conds = sprintf("%s and filename not like  \'\%%%s\%\'",$conds,$item);
    }
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
    if ($tp eq "G4Hits" || $tp eq "G4HitsOld")
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
	    my $splitstring = sprintf("_%s",$pileupstring);
            my @sp2 = split(/$splitstring/,$filenamestring_with_runnumber);
	    $newfilenamestring = sprintf("%s-%010d-",$sp2[0],$runnumber);
	}
	my $newgetfilesql = $getfilesql;
	$newgetfilesql =~ s/$filenamestring_with_runnumber/$newfilenamestring/;
	$getfiles{"G4Hits"} = $dbh->prepare($newgetfilesql);
	$getfiles{"G4HitsOld"} = $dbh->prepare($newgetfilesql);
	if (defined $verbose)
	{
	    print "sql (newgetfilesql): $newgetfilesql\n";
	}
    }
    else
    {
	$getfiles{$tp} = $dbh->prepare($getfilesql);
	if (defined $verbose)
	{
	    print "sql (getfilesql): $getfilesql\n";
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

my $entries = 200000000; # given that we have 200k files max, this value is always higher
my $lowtype;
# here we find the dst type with the smallest number of entries (segments)
# so we do not loop too much when finding matches for the other types
if (defined $verbose)
{
    print "hashing done, finding hash with lowest number of entries\n";
}
foreach my $tp (sort { $a <=> $b } keys %allfilehash)
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
foreach my $seg (sort { $a <=> $b } keys %{$allfilehash{$lowtype}})
{
    foreach my $tp (sort { $a <=> $b } keys %allfilehash)
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
my %filesorted = ();
foreach my $seg (@segarray)
#foreach my $seg (sort { $a <=> $b } @segarray)
{
    $nSelectedEvents += $allevthash{$lowtype}{$allfilehash{$lowtype}{$seg}};
#	print "segment $seg is good\n";
#    foreach my $tp (sort keys %allfilehash)
    foreach my $tp (keys %allfilehash)
    {
#	    print "using $allfilehash{$tp}{$seg}\n";
#	my $printcmd = sprintf("echo %s >> %s",$allfilehash{$tp}{$seg},$dsttype{$tp});
	$filesorted{$dsttype{$tp}}{$allfilehash{$tp}{$seg}} = 1;
#	print "$printcmd\n";
#	system($printcmd);
    }

}
foreach my $listfile (keys %filesorted)
{
    open(F3,">$listfile");
    foreach my $fil (sort keys %{$filesorted{$listfile}})
    {
	print F3 "$fil\n";
    }
    close(F3);
}
print "wrote the following list files containing >= $nSelectedEvents events:\n";
foreach my $tp (sort { $a <=> $b } keys %allfilehash)
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
# for no pileup pass(X) --> pass(X-1)
# pass1
    $filetypes{"G4Hits"} = "G4 Hits";
#    $filetypes{"G4HitsOld"} = "Old G4 Hits";
# pass2
    $filetypes{"DST_BBC_G4HIT"} = "Pileup BBC (now MBD), EPD G4Hits";
    $filetypes{"DST_CALO_G4HIT"} = "Pileup Calorimeter G4Hits";
    $filetypes{"DST_TRKR_G4HIT"} = "Pileup Tracking Detector G4 Hits";
    $filetypes{"DST_TRUTH_G4HIT"} = "temporary Pileup Truth info, use DST_TRUTH";
# pass3 mbdepd
    $filetypes{"DST_MBD_EPD"} = "Reconstructed Mbd, Epd";
# pass3 calo
    $filetypes{"DST_CALO_CLUSTER"} = "Reconstructed Calorimeter Towers and Clusters";
#pass3 trk
    $filetypes{"DST_TRKR_HIT"} = "TPC and Silicon Hits";
    $filetypes{"DST_TRUTH"} = "Truth Info (updated with Clusters)";
#pass4 truth jets
    $filetypes{"DST_TRUTH_JET"} = "Truth Jets";
#pass4 tracks
    $filetypes{"DST_TRKR_CLUSTER"} = "pass0 output: tpc clusters";
    $filetypes{"DST_TRACKSEEDS"} = "passA output: track seeds";
    $filetypes{"DST_TRACKS"} = "passC output: Reconstructed Tracks";
#pass5 tracks/clusters
    $filetypes{"DST_GLOBAL"} = "Global Info (MBD, sEPD, Vertex)";
    $filetypes{"DST_TRUTH_RECO"} = "digested track truth info";
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
	if ($sp1[1] =~ /$runsplit_embed/)
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
	    my $part = $1;
            my $mom;
            my $minp = $2;
            my $maxp = $3;
	    if ($part =~ /(\S+)_(\S+)/)
	    {
		$part = $1;
		$mom = $2;
	    }
	    if (defined $mom)
	    {
		print "CreateFileList.pl -type 14 $types{$name} -run $runnumber -particle $part -mom $mom -pmin $minp -pmax $maxp\n";
	    }
	    else
            {
                print "CreateFileList.pl -type 14 $types{$name} -run $runnumber -particle $part -pmin $minp -pmax $maxp\n";
            }
	}
        else
        {
            print "CreateFileList.pl -type 14 $types{$name} -run $runnumber -particle $name\n";

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

sub print_runs
{
    my $getrunnumbers = $dbh->prepare("select distinct(runnumber) from datasets where dataset = 'mdc2' order by runnumber");
    $getrunnumbers->execute();
    print "Available Runs (check our wiki for more details for each runnumber):\n";
    while(my @res = $getrunnumbers->fetchrow_array())
    {
	print "$res[0]\n";
    }
    print "NB: Not all DSTs are available for all runs\n";
    $getrunnumbers->finish();
}
