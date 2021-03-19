#!/usr/bin/perl
use lib './';
use Getopt::Std;
use Sys::Hostname;

getopts("cde:f:g:hmn:oq:rs:tvwz");

# Print help message if -h option is specified
if($opt_h) {
    print "Usage: pov-movie.pl [-options <n> ...] SIM-DIR\n\n";
    print "Options:\n";
    print "-c             (Add cylinders)\n";
    print "-d             (Don't duplicate existing files)\n";
    print "-e <n>         (Only render every <num> frame)\n";
    print "-f <n>         (Start redener from <num> frame)\n";
    print "-g <n>         (Render up to <num> frame)\n";
    print "-h             (Print this information)\n";
    print "-m             (Switch off the mesh)\n";
    print "-o             (Render mesh as triangles, without normals)\n";
    print "-q <n>         (Quality of rendering, 1=good, 3=extreme)\n";
    print "-s <n>         (Render a single frame)\n";
    print "-t             (Add tracers)\n";
    print "-v             (Verbose output)\n";
    print "-w             (Disable making a movie)\n";
    print "Example:\nperl pov-movie.pl -t -d -q 3 mySim.out\n";
    print "Renders any unrendered output (option -d) files in mySim.out in the highest quality (option -q 3), with tracers (option -t).\n";
    print "Tracers, as well as meshes of the solid objects are rendered.\n";
    print "Mesh files (ctr.%05d) and tracer files (tr.%05d) must exsit in the sim directory.\n";

    exit 0;
}

die "One or two arguments required" unless @ARGV==1;

# Set miscellaneous variables, and those used to control which frames are
# rendered
$dr=$ARGV[0];
$msh=$opt_o?"mtr":"msh";
$verb=$opt_v?"":">/dev/null 2>/dev/null";
$every=$opt_e?$opt_e:1;
$g=defined $opt_g?$opt_g:-1;

# Read the first line of the POV-Ray header file to get the rendering
# dimensions. Assemble the POV-Ray flags.
$opt_q=1 unless defined $opt_q;
die "POV quality out of range\n" if $opt_q<0 || $opt_q>3;
open B,"../pov_headers/production.pov" or die "Can't open POV-Ray header file\n";
$_=<B>;
m/^\/\/ W=(\d*) H=(\d*)/ or die "Can't read rendering size\n";
$pov_opts="+W$1 +H$2 ".
          ("","+R3 +A0.01 -J","+R6 +A0.001 -J","+R9 +A0.0001 -J")[$opt_q];

# Read, process, and store the rest of the POV-Ray header file
$gpn=0;
while(<B>) {
    if($opt_c) {s/^CYL://;} else {next if /^CYL:/;}
    if($opt_m) {next if /^MSH:/;} else {s/^MSH://;}
    if($opt_t) {s/^TRA://;} else {next if /^TRA:/;}
    $gp[$gpn++]=$_;
}
close B;$gpn--;

# Loop over the available frames
$a=defined $opt_s?$opt_s:(defined $opt_f?$opt_f : 0);
$ctr_fn=sprintf("%s/ctr.%05d", $dr, $a);
while(-e $ctr_fn and ($g==-1 or $a<=$opt_g)) {

    # Assemble output filename and check for skipping/termination conditions
    $fn=sprintf "fr_%05d.png",$a;
    last if defined $opt_s && $a>$opt_s;
    do {
	    #print("Skip CTR $ctr_fn \n");
	    #print("Skip A $a\n");
        $a++;
        $ctr_fn=sprintf("%s/ctr.%05d", $dr, $a);
	next;
    } if defined $opt_d && -e "$dr/$fn";

    # Assemble the POV file
    $ctr_fn=sprintf("%s/ctr.%05d", $dr, $a);
    #print("CTR $ctr_fn \n");
    #print("A $a\n");
    $tr_fn=sprintf("%s/tr.%05d", $dr, $a);
    $tf="rtemp.pov";
    open A,$opt_r?"| bzip2 -9 -c >$dr/$tf.bz2":">$dr/$tf" or die "Can't open temporary POV file\n";
    foreach $i (0..$gpn) {
        $_=$gp[$i];

        if(/^#include "msh\.pov"/) {
            print A `./unpack $msh $ctr_fn -`;
        } elsif(/^#include "cyl\.pov"/) {
            print A `./unpack cyl $ctr_fn -`;
        } elsif(/^#include "sph\.pov"/) {
            print A `./tr_unpack sph $tr_fn -`;
        } else {
            print A;
        }
    }
    close A;

    # Render the frame, forking jobs to remote processors if the "-r" option is
    # given

    $pov_cmd="nice -n 19 povray $tf -D +O$fn $pov_opts";

    # Run POV-Ray locally
    print "Frame $a\n";
    system "cd $dr; $pov_cmd $verb";

    $a+=$every;
    $ctr_fn=sprintf("%s/ctr.%05d", $dr, $a);
}

# Additional code to automatically make a movie
unless ($opt_w) {
    ($mf=$dr)=~s/\.out//;
    unlink "$mf.mov";
    system "ffmpeg -r 30 -y -i $dr/fr_%05d.png -preset veryslow -c:v libx265 -crf 17 -pix_fmt yuv420p -tag:v hvc1 -movflags faststart $mf.mov";
}
