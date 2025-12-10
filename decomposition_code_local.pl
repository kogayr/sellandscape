#!/usr/bin/perl
############################################################
# Template: 1.5, August 09, 2018
# Parent Path: /home/wolf/bin on frosty
############################################################

############################################################
# Replacements for YIW modules
############################################################
{
    package YIW::basic;
    use strict;
    use warnings;
    our (@myGlobList, %myOptList);

    sub my_args {
        my ($rargs, $need) = @_;
        $need //= 0;
        @myGlobList = ();
        %myOptList  = ();
        for (my $i = 0; $i < @$rargs; $i++) {
            my $x = $$rargs[$i];
            if ($x =~ /^-(\w+)(=(.*))?$/) {
                my ($k, $v) = ($1, defined $3 ? $3 : 1);
                $myOptList{$k} = $v;
            } else {
                push @myGlobList, $x;
            }
        }
        die "Usage error: missing argument\n" if ($need && !@myGlobList);
        return 1;
    }
}

{
    package YIW::stat;
    use strict;
    use warnings;

    sub match_norm_miq {
        my ($arrref) = @_;
        my @vals = ref($arrref) eq 'ARRAY' ? @{$arrref} : ();
        die "match_norm_miq: empty input" unless @vals;
        my @s = sort { $a <=> $b } @vals;
        my $n = scalar @s;
        my $median = ($n % 2) ? $s[int($n/2)] : ($s[$n/2 - 1] + $s[$n/2]) / 2;
        my @absd = sort { $a <=> $b } map { abs($_ - $median) } @s;
        my $mad = ($n % 2) ? $absd[int($n/2)] : ($absd[$n/2 - 1] + $absd[$n/2]) / 2;
        my $sigma = 1.482602218505602 * $mad;
        if ($sigma == 0) {
            my ($sum, $sumsq) = (0, 0);
            foreach my $v (@s) { $sum += $v; $sumsq += $v*$v; }
            my $mean = $sum / $n;
            my $var  = ($sumsq / $n) - ($mean*$mean);
            $var = 0 if $var < 0 && $var > -1e-15;
            $sigma = sqrt($var);
        }
        $sigma = 1e-300 if $sigma == 0;
        return ($median, $sigma);
    }

    sub max { $_[0] > $_[1] ? $_[0] : $_[1] }
}

{
    package main;
    use Math::CDF qw(pnorm);

    sub int_commify {
        local $_ = shift;
        1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
        return $_;
    }

    our (@myGlobList, %myOptList);
    *myGlobList = \@YIW::basic::myGlobList;
    *myOptList  = \%YIW::basic::myOptList;
    *max        = \&YIW::stat::max;
}

############################################################
# System etc
############################################################

# uncomment the next line if better randomness is important
# use Math::Random::MT qw/rand srand/;

our $Pver = "0.1";
our $Pdat = "2024-09-01";
our ($Pnam) = ($0 =~ m/([^\/]+)$/);
our ($Path) = ($0 =~ m/^(.+)\/[^\/]+$/);
$Path = "." unless ($Path);
our $Ppid = $$; # Process ID
our $Base = "tmp.$Pnam.$Ppid";

{
    my ($host) = ($ENV{"HOSTNAME"} =~ m/^(\w+)/);
    my $salt = sprintf ".%s.%03d", $host, int(rand(1000));
    $Base .= $salt;
}

our $CMD = "";
our $DEBUG = 0;
our $VERBOSE = 1;
our $DOIT = 1;
our $EXITCODE = 0;

############################################################
# Definitions
############################################################
my $EPSILON_e2 = 1e-6;
my $ITER_max   = 1000;
my $EVALUE_thr = 0.5;
my $cleanup    = 0;

############################################################
# Global variables
############################################################

############################################################
# Instructions etc
############################################################
$Instructions = <<EOINPUT;
$Path/$Pnam $Pver, $Pdat

Use:
    $Pnam name [options]

Options:
    -clean  dismiss points with unlikely residuals and redo
EOINPUT

############################################################
# code start
############################################################

#--- get and process arguments --------------------------
YIW::basic::my_args(\@ARGV, 1);  # pass 0 to block STDIN
!@myGlobList and print $Instructions and exit 0;

if (exists $myOptList{"DEBUG"}) {   # hidden common option
    $DEBUG   = 1;
    $VERBOSE = 0x7fff;
}
$VERBOSE = int $myOptList{"VERB"} if (exists $myOptList{"VERB"});
$DOIT = 0 if (exists $myOptList{"IDLE"});
$Base = $myOptList{"BASE"} if ($myOptList{"BASE"} ne "");
$cleanup = 1 if (exists $myOptList{"clean"});

my $fdat = shift @myGlobList;

#--- rest of the code -----------------------------------
if ($CMD ne "") {
    print STDERR "$Pnam:\t[ $CMD ]\n" if ($VERBOSE);
    $EXITCODE = system $CMD if ($DOIT);
    die "$Pnam:\t[ $CMD ] failed with [ $EXITCODE ]" if ($EXITCODE != 0);
}

$; = "\t";
my %dnds  = ();
my @latgc = ();
my @lcogs = ();

{
    my %tmpa = ();
    my %tmpc = ();
    print STDERR "$Pnam:\tReading $fdat\n" if ($VERBOSE);
    open HAND, "<$fdat" or die "$Pnam:\tCan't read \"$fdat\"";
    while (<HAND>) {
        chomp;
        my ($atgc, $ccog, $xx) = split /\t/;
        $dnds{$atgc, $ccog} = log($xx);
        $tmpa{$atgc}++;
        $tmpc{$ccog}++;
    }
    close HAND;
    @latgc = sort keys %tmpa;
    @lcogs = sort keys %tmpc;
}

printf STDERR "$Pnam:\t\t%s\tATGCs\n", int_commify(scalar @latgc) if ($VERBOSE);
printf STDERR "$Pnam:\t\t%s\tCOGs\n", int_commify(scalar @lcogs) if ($VERBOSE);
printf STDERR "$Pnam:\t\t%s\tdatapoints\n", int_commify(scalar keys %dnds) if ($VERBOSE);

#--- decompose 1st pass ---------------------------------
my @qatgc = ();
my @qcogs = ();
my %ddres = ();
decompose_dataset(\%dnds, \@latgc, \@lcogs, \@qatgc, \@qcogs, \%ddres);

#--- print 1st pass -------------------------------------
if ($cleanup <= 0) {
    printf "# ATGC\n";
    for (my $i = 0; $i < @latgc; $i++) {
        printf "%s\t%f\n", $latgc[$i], $qatgc[$i];
    }
    printf "# COGs\n";
    for (my $i = 0; $i < @lcogs; $i++) {
        printf "%s\t%f\n", $lcogs[$i], $qcogs[$i];
    }
    printf "# residuals\n";
    foreach my $xx (sort keys %ddres) {
        printf "%s\t%f\n", $xx, $ddres{$xx};
    }
    exit;
}

#--- residuals stats ------------------------------------
my @lres = values %ddres;
my ($noav, $nosd) = YIW::stat::match_norm_miq(\@lres);
printf STDERR "$Pnam:\t\t%f\t%f\tresiduals normal equivalent\n", $noav, $nosd if ($VERBOSE);

#--- reassess -------------------------------------------
my %dxdx  = ();
my @latgx = ();
my @lcogx = ();
my $ndat  = scalar keys %dnds;
my $ndis  = 0;

{
    my %tmpa = ();
    my %tmpc = ();
    foreach my $xx (keys %dnds) {
        my $zz = ($ddres{$xx} - $noav) / $nosd;
        my $ev = pnorm(-abs($zz)) * $ndat;  # <-- replaced cdf_norm_s() with pnorm()
        if ($ev < $EVALUE_thr) {
            $ndis++;
        } else {
            $dxdx{$xx} = $dnds{$xx};
            my ($atgc, $ccog) = split /$;/, $xx;
            $tmpa{$atgc}++;
            $tmpc{$ccog}++;
        }
    }
    @latgx = sort keys %tmpa;
    @lcogx = sort keys %tmpc;
}

printf STDERR "$Pnam:\t\t%s\tdismissed\n", int_commify($ndis) if ($VERBOSE);
printf STDERR "$Pnam:\t\t%s\tATGCs\n", int_commify(scalar @latgx) if ($VERBOSE);
printf STDERR "$Pnam:\t\t%s\tCOGs\n", int_commify(scalar @lcogx) if ($VERBOSE);
printf STDERR "$Pnam:\t\t%s\tdatapoints\n", int_commify(scalar keys %dxdx) if ($VERBOSE);

#--- decompose 2nd pass ---------------------------------
my @qatgx = ();
my @qcogx = ();
my %ddrex = ();
decompose_dataset(\%dxdx, \@latgx, \@lcogx, \@qatgx, \@qcogx, \%ddrex);

#--- print 2nd pass -------------------------------------
if ($cleanup > 0) {
    printf "# ATGC\n";
    for (my $i = 0; $i < @latgx; $i++) {
        printf "%s\t%f\n", $latgx[$i], $qatgx[$i];
    }
    printf "# COGs\n";
    for (my $i = 0; $i < @lcogx; $i++) {
        printf "%s\t%f\n", $lcogx[$i], $qcogx[$i];
    }
    printf "# residuals\n";
    foreach my $xx (sort keys %ddrex) {
        printf "%s\t%f\n", $xx, $ddrex{$xx};
    }
    exit;
}

#--- clean ----------------------------------------------
# unlink <$Base.*> unless ($DEBUG);

############################################################
# decompose_dataset($rdat,$rdi1,$rdi2,$qdi1,$qdi2,$rres)
############################################################
sub decompose_dataset {
    my ($rdat, $rdi1, $rdi2, $qdi1, $qdi2, $rres) = @_;
    @$qdi1 = 0 x @$rdi1;
    @$qdi2 = 0 x @$rdi2;
    %$rres = %$rdat;

    my $e0 = calc_sumsq($rres);
    my $e1 = $e0;
    my $dd = 1;
    printf STDERR "$Pnam:\t\tIt %d\t%f\t%f\n", 0, $dd, $e1 if ($VERBOSE);

    my $iter = 1;
    while (1) {
        my @ave2 = ();
        calc_ave2($rres, $rdi1, $rdi2, \@ave2);
        for (my $j = 0; $j < @$rdi2; $j++) {
            $$qdi2[$j] += $ave2[$j];
        }
        calc_residuals($rdat, $rdi1, $rdi2, $qdi1, $qdi2, $rres);

        my @ave1 = ();
        calc_ave1($rres, $rdi1, $rdi2, \@ave1);
        for (my $i = 0; $i < @$rdi1; $i++) {
            $$qdi1[$i] += $ave1[$i];
        }
        for (my $i = 0; $i < @$rdi1; $i++) {
            $$qdi1[$i] -= $$qdi1[0];
        }
        calc_residuals($rdat, $rdi1, $rdi2, $qdi1, $qdi2, $rres);

        $e1 = calc_sumsq($rres);
        $dd = abs($e1 - $e0) / max($e1, $e0);
        if ($iter % 50 == 0) {
            printf STDERR "$Pnam:\t\tIt %d\t%f\t%f\n", $iter, $dd, $e1 if ($VERBOSE);
        }
        last if ($dd <= $EPSILON_e2);
        $e0 = $e1;
        $iter++;
        last if ($iter > $ITER_max);
    }
    printf STDERR "$Pnam:\t\tIt %d\t%f\t%f\n", $iter, $dd, $e1 if ($VERBOSE);
}

############################################################
# calc_sumsq($rdat)
############################################################
sub calc_sumsq {
    my $rdat = shift;
    my $s2 = 0;
    foreach my $yy (values %$rdat) {
        $s2 += $yy * $yy;
    }
    return $s2;
}

############################################################
# calc_ave1($rdat,$rdi1,$rdi2,$rav1)
############################################################
sub calc_ave1 {
    my ($rdat, $rdi1, $rdi2, $rav1) = @_;
    for (my $i = 0; $i < @$rdi1; $i++) {
        my $dim1 = $$rdi1[$i];
        my $ss = 0;
        my $nn = 0;
        for (my $j = 0; $j < @$rdi2; $j++) {
            my $dim2 = $$rdi2[$j];
            next unless (exists $$rdat{$dim1, $dim2});
            $ss += $$rdat{$dim1, $dim2};
            $nn++;
        }
        $ss /= $nn if ($nn > 0);
        $$rav1[$i] = $ss;
    }
}

############################################################
# calc_ave2($rdat,$rdi1,$rdi2,$rav2)
############################################################
sub calc_ave2 {
    my ($rdat, $rdi1, $rdi2, $rav2) = @_;
    for (my $i = 0; $i < @$rdi2; $i++) {
        my $dim2 = $$rdi2[$i];
        my $ss = 0;
        my $nn = 0;
        for (my $j = 0; $j < @$rdi1; $j++) {
            my $dim1 = $$rdi1[$j];
            next unless (exists $$rdat{$dim1, $dim2});
            $ss += $$rdat{$dim1, $dim2};
            $nn++;
        }
        $ss /= $nn if ($nn > 0);
        $$rav2[$i] = $ss;
    }
}

############################################################
# calc_residuals($rdat,$rdi1,$rdi2,$qdi1,$qdi2,$rres)
############################################################
sub calc_residuals {
    my ($rdat, $rdi1, $rdi2, $qdi1, $qdi2, $rres) = @_;
    for (my $i = 0; $i < @$rdi1; $i++) {
        my $dim1 = $$rdi1[$i];
        for (my $j = 0; $j < @$rdi2; $j++) {
            my $dim2 = $$rdi2[$j];
            next unless (exists $$rdat{$dim1, $dim2});
            $$rres{$dim1, $dim2} = $$rdat{$dim1, $dim2} - $$qdi1[$i] - $$qdi2[$j];
        }
    }
}
