#!/usr/local/bin/perl

########################################################################
# Copyright 1998-2020 CERN for the benefit of the EvtGen authors       #
#                                                                      #
# This file is part of EvtGen.                                         #
#                                                                      #
# EvtGen is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# EvtGen is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     #
########################################################################

# B0 ->  K+ pi- pi0
# B+ -> pi+ pi+ pi-
# B+ -> pi+ pi0 pi0
# B0 -> pi+ pi- pi0
#
# Reference: A. Snyder and H. Quinn [Phys.Rev.D48 (1993) 2139]
#
# This model was implemented in dedicated FORTRAN executables. 
# The parameters used are the ones resulting from the ALEPH analysis
# which might be found in CERN-PPE:97013 
#
# This script generates decay files for a generic c++ model PTO3P(_CP).
# Strong phase is set to zero and not considered.

use Getopt::Long;
use Math::Complex;


# Subroutines
#____________

sub print_res {

    my $rH    = shift;
    my $conj  = shift;
    my $res  = $$rH{"res"};
    my $part = $$rH{"part"};
    my $m0   = $$rH{"m0"};
    my $g0   = $$rH{"g0"};
    my $c    = $$rH{"c"};
    my $ang  = $$rH{"ang"};

    my $re = Re($c);
    my $im = Im($c);

    $im = -$im if ($conj eq "yes");

    print  "RESONANCE\t$res\t$part\t$m0 $g0\n";
    printf "CARTESIAN\t%4.4f %4.4f\t$ang\n", $re, $im;
    print  "TYPE\t\tRBW_KUEHN\n\n";
}

sub print_res_array {

    my $rA = shift;
    my $conj = shift;
    for my $i (0..$#$rA) {
	
	my $rH = $$rA[$i];
	print_res($rH,$conj);
    }
}

sub print_decay {

    my $rH = shift;
    my $mother = $$rH{"mother"}; 
    my $daus = $$rH{"daus"};
    my $model = $$rH{"model"}; 
    my $maxpdf = $$rH{"maxpdf"}; 
    my $rA = $$rH{"res"};
    print "Decay $mother\n\n";
    print "1.0 $daus $model\n\n";

    if($model eq "PTO3P_CP") {
	print "dm \n\n";
    }

    print "MAXPDF $maxpdf\n\n";
	
    print_res_array($rA,"no");

    if($model eq "PTO3P_CP") {
	print "CONJUGATE\n\n";
	print_res_array($rA,"yes");
    }

    print ";\n";
    print "Enddecay\n";
    print "End\n";
}


# CKM matrix parameters and parameters of the Snyder-Quinn model 
#_______________________________________________________________

$alpha = 1.365;
$beta  = 0.39;
$ema = cos($alpha) - i*sin($alpha);
$epb = cos($beta)  + i*sin($beta);
$rho = 0.;
$eta = 0.5;
$c_f = sqrt(($rho*$rho+$eta*$eta)/((1.-$rho)*(1.-$rho)+$eta*$eta));

$tp0 = 1.19 * $ema * $c_f;
$tm0 = 1.19 * $ema * $c_f;
$t0m = 0.8  * $ema * $c_f;
$t0p = 0.8  * $ema * $c_f;
$tpm = 0.59 * $ema * $c_f;
$tmp = 1.08 * $ema * $c_f;
$p1  = -0.02;
$p0  = 0.03;

# B->3pi matrix elements
$m1  = $tp0 + 2*$p1;
$m2  = $t0p - 2*$p1;
$m3  = $tpm + $p1 + $p0;
$m4  = $tmp - $p1 + $p0;
$m5  = -$tpm - $tmp + $tp0 + $t0p - 2.*$p0;
$n1 = ~$m1;
$n2 = ~$m2;
$n3 = ~$m3;
$n4 = ~$m4;
$n5 = ~$m5;

# rho, rho', rho'' Fractions
$b = -0.229;
$g =  0.075;

# B0->Kpipi matrix elements
$m_kstm = 0.22*$ema - 1.2*$epb;
$m_kst0 = 0.015*$ema + 0.850*$epb;
$m_rho  = 0.13*$ema + 0.16*$epb;

# masses
$mrho = 0.7734; $grho = 0.1477;
$mrhop = 1.465; $grhop = 0.696;
$mrhopp = 1.760; $grhopp = 0.215;


# Data structure for resonances
#______________________________

%r11 = ( res => "AC", part => "K*-", m0 => 0.89, g0 => 0.0498, c => $m_kstm, ang => "AB" );
%r12 = ( res => "AB", part => "anti-K*0", m0 => 0.8961, g0 => 0.0505, c => $m_kst0, ang => "AC" );
%r13 = ( res => "BC", part => "rho+", m0 => 0.770, g0 => 0.150, c => $m_rho, ang => "AB" );
@r1  = ( \%r11, \%r12, \%r13 );
%h1 = ( mother => "B0", daus =>  "K+ pi- pi0", model => "PTO3P", res => \@r1, maxpdf => 520. );

%r211 = ( res => "AC", part => "rho0", m0 => $mrho, g0 => $grho, c => $m2/(1.+$b+$g), ang => "BC" );
%r212 = ( res => "AC", part => "rho(2S)0", m0 => $mrhop, g0 => $grhop, c => $m2*$b/(1.+$b+$g), ang => "BC" );
%r213 = ( res => "AC", part => "rho(3S)0", m0 => $mrhopp, g0 => $grhopp, c => $m2*$g/(1.+$b+$g), ang => "BC" );
%r221 = ( res => "BC", part => "rho0", m0 => $mrho, g0 => $grho, c => $m2/(1.+$b+$g), ang => "AC" );
%r222 = ( res => "BC", part => "rho(2S)0", m0 => $mrhop, g0 => $grhop, c => $m2*$b/(1.+$b+$g), ang => "AC" );
%r223 = ( res => "BC", part => "rho(3S)0", m0 => $mrhopp, g0 => $grhopp, c => $m2*$g/(1.+$b+$g), ang => "AC" );
@r2  = ( \%r211, \%r212, \%r213, \%r221, \%r222, \%r223 );
%h2 = ( mother => "B+", daus => "pi+ pi+ pi-", model => "PTO3P", res => \@r2, maxpdf => 7.5 );

%r311 = ( res => "AB", part => "rho+", m0 => $mrho, g0 => $grho, c => $m1/(1.+$b+$g), ang => "AC" );
%r312 = ( res => "AB", part => "rho(2S)+", m0 => $mrhop, g0 => $grhop, c => $m1*$b/(1.+$b+$g), ang => "AC" );
%r313 = ( res => "AB", part => "rho(3S)+", m0 => $mrhopp, g0 => $grhopp, c => $m1*$g/(1.+$b+$g), ang => "AC" );
%r321 = ( res => "AC", part => "rho+", m0 => $mrho, g0 => $grho, c => $m1/(1.+$b+$g), ang => "AB" );
%r322 = ( res => "AC", part => "rho(2S)+", m0 => $mrhop, g0 => $grhop, c => $m1*$b/(1.+$b+$g), ang => "AB" );
%r323 = ( res => "AC", part => "rho(3S)+", m0 => $mrhopp, g0 => $grhopp, c => $m1*$g/(1.+$b+$g), ang => "AB" );
@r3  = ( \%r311, \%r312, \%r313, \%r321, \%r322, \%r323 );
%h3 = ( mother => "B+", daus => "pi+ pi0 pi0", model => "PTO3P", res => \@r3, maxpdf => 15.0 );

%r411 = ( res => "AC", part => "rho+", m0 => $mrho, g0 => $grho, c => $m3/(1.+$b+$g), ang => "AB" );
%r412 = ( res => "AC", part => "rho(2S)+", m0 => $mrhop, g0 => $grhop, c => $m3*$b/(1.+$b+$g), ang => "AB" );
%r413 = ( res => "AC", part => "rho(3S)+", m0 => $mrhopp, g0 => $grhopp, c => $m3*$g/(1.+$b+$g), ang => "AB" );
%r421 = ( res => "BC", part => "rho-", m0 => $mrho, g0 => $grho, c => $m4/(1.+$b+$g), ang => "AC" );
%r422 = ( res => "BC", part => "rho(2S)-", m0 => $mrhop, g0 => $grhop, c => $m4*$b/(1.+$b+$g), ang => "AC" );
%r423 = ( res => "BC", part => "rho(3S)-", m0 => $mrhopp, g0 => $grhopp, c => $m4*$g/(1.+$b+$g), ang => "AC" );
%r431 = ( res => "AB", part => "rho0", m0 => $mrho, g0 => $grho, c => $m5/(1.+$b+$g), ang => "AC" );
%r432 = ( res => "AB", part => "rho(2S)0", m0 => $mrhop, g0 => $grhop, c => $m5*$b/(1.+$b+$g), ang => "AC" );
%r433 = ( res => "AB", part => "rho(3S)0", m0 => $mrhopp, g0 => $grhopp, c => $m5*$g/(1.+$b+$g), ang => "AC" );

@r4  = ( \%r411, \%r412, \%r413, \%r421, \%r422, \%r423, \%r431, \%r432, \%r433);
%h4 = ( mother => "B0", daus => "pi+ pi- pi0", model => "PTO3P_CP", res => \@r4, maxpdf => 17.9 );


# Get command line options and 
# print a decay file

my $i;
my $result = GetOptions('n=i' => \$i);

   if($i == 1) { print_decay(\%h1); break; }
elsif($i == 2) { print_decay(\%h2); break; }
elsif($i == 3) { print_decay(\%h3); break; }
elsif($i == 4) { print_decay(\%h4); break; }
else           { print "Usage: evt_bdalitz.pl -n <Decay Number>\n"; }


