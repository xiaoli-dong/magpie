
package Constants;

$revcomp = -4;
$forward = 4; 
$nostrand = 0;

my $base = $forward;
$dna = ++$base;
$protein = ++$base;
$motif = ++$base;
$trna  = ++$base;
$model = ++$base;
$hmm = ++$base;
$msa = ++$base;
$repeat= ++$base;
$orf = ++$base;

$localtool = ++$base;
$httptool = ++$base;
$emailtool = ++$base;
$facttool = ++$base;
$texttool = ++$base;

$allscope = ++$base;
$dnascope = ++$base;
$aascope = ++$base;
$orfscope = ++$base;
$interorfscope = ++$base;
$groupscope = ++$base;
$projectscope = ++$base;

$private = ++$base;
$public = ++$base;

$dbid = ++$base;
$word = ++$base;
$seqid = ++$base;  # A Magpie identifier
$goid = ++$base;  # Gene Ontology ID
1;
__END__

