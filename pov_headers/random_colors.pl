#!/usr/bin/perl

$fn=$ARGV[0];
open A,$fn or die "Can't open input file\n";

$fn=~s/\.in\.pov$/\.pov/g or die "Must run on '.in.pov' file\n";
open B,">$fn" or die "Can't open output file\n";

while(<A>) {
	if(/^\/\/ randcols/) {
		(undef,undef,$s1,$s2,$hs,$he,$ss,$se,$vs,$ve)=split;
		foreach $i ($s1..$s2) {
			$h=$hs+($he-$hs)*rand;
			$s=$ss+($se-$ss)*rand;
			$v=$vs+($ve-$vs)*rand;

			# Ensure hue between 0 and 1
			$hi=int $h;
			$h-=$hi;

			$c=$v*$s;
			$h*=6;
			if($h<1) {
				$r=$c;$g=$h*$c;$b=0;
			} elsif($h<2) {
				$r=(2-$h)*$c;$g=$c;$b=0;
			} elsif($h<3) {
				$r=0;$g=$c;$b=($h-2)*$c;
			} elsif($h<4) {
				$r=0;$g=(4-$h)*$c;$b=$c;
			} elsif($h<5) {
				$r=($h-4)*$c;$g=0;$b=$c;
			} else {
				$r=$c;$g=0;$b=(6-$h)*$c;
			}
			$r+=$v-$c;$g+=$v-$c;$b+=$v-$c;
			printf B "#declare t_msh$i=texture{pigment{rgb <%.3f,%.3f,%.3f>} finish{f0}}\n",$r,$g,$b;
		}
	} else {
		print B;
	}
}
