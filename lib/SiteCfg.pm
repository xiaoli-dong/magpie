
package SiteCfg;

my $fastadir = $ENV{MAGPIEHOME}."/ext-bin";
my $blastdir = "/usr/local/blast";
my $blastdb  = "/export/data/databases/blast/";
$estscandir = "/export/data/programs/ESTScan3";
$readseq  = $ENV{MAGPIEHOME}."/ext-bin/readseq";
$extract_seq  = $ENV{MAGPIEHOME}."/ext-bin/extract_seq";
$gzip = "/bin/gzip";
$tar  = "/bin/tar";
$ps   = "/bin/ps";

$fasta = "$fastadir/fasta34_t";
$fastx = "$fastadir/fastx34_t";
$fasts = "$fastadir/fasts34_t";
$fastf = "$fastadir/fastf34_t";
$fasty = "$fastadir/fasty34_t";
$tfasty = "$fastadir/tfasty34_t";
$tfastf = "$fastadir/tfastf34_t";
$tfastx = "$fastadir/tfastx34_t";
$tfasts = "$fastadir/tfasts34_t";

$megablast = "$blastdir/megablast";
$blastall = "$blastdir/blastall";
$rpsblast = "$blastdir/rpsblast";
$psiblast = "$blastdir/psiblast";
$formatdb = "$blastdir/formatdb";
$blastmat = "$blastdir/data";

1;

__END__
