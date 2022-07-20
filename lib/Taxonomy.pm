
package Taxonomy;

use GDBM_File;
use FileHandle;
use Fcntl;
use LWP::Simple;
use vars qw(%hierarchy);

$login = getlogin();

#Top level taxids
$viruses  = 10239;
$bacteria = 2;
$archaea  = 2157;
$eukarya  = 2759;
$plants   = 33090;
$metazoa  = 33208;
$fungi    = 4751;
$metagenomes = 408169;
my $i;
%hierarchy = (
superkingdom => $i++,
kingdom => $i++,
superphylum => $i++,
phylum => $i++,
subphylum => $i++,
superclass => $i++,
infraclass => $i++,
class => $i++,
subclass => $i++,
superorder => $i++,
order => $i++,
suborder => $i++,
infraorder => $i++,
parvorder => $i++,
superfamily => $i++,
family => $i++,
subfamily => $i++,
tribe => $i++,
subtribe => $i++,
genus => $i++,
subgenus => $i++,
"species group" => $i++,
"species subgroup" => $i++,
species => $i++,
subspecies => $i++,
forma => $i++,
varietas => $i++,
"no rank" => $i++
);

#Set this to the location of the NDBM index files for the taxonomy database
my $species_index_path = $ENV{'MAGPIEHOME'}."/Tables/Taxonomy/index";

sub new{
  return bless {};
}

#The lower the number, the higher in the taxonomy hierarchy this node is
sub rank_order{
  return $hierarchy{$_[1]};
}

sub lineage{
  my ($self, @lineage) = @_;
  my $node = $lineage[0];


  while($node = $self->parent($node)){


    last if $node == $lineage[$#lineage] || !defined $node; #Circular reference (e.g. no rank)
    push @lineage, $node;
  }
  return reverse @lineage;
}

# Calculates a "distance" between taxonomic identifiers as an approximation of evolutionary distance
sub dist{
  my ($self, $taxid_a, $taxid_b) = @_;

  my @lineage_a = $self->lineage($taxid_a);
  my @lineage_b = $self->lineage($taxid_b);
  my $dist = 0;
  my $i = 0;
  SEARCH:{for my $tax_a (reverse @lineage_a){
    for my $tax_b (reverse @lineage_b){
      if($tax_b == $tax_a){
        $dist = $i;
        last SEARCH;
      }
    }
    $i++;
  }
  }
  return $dist;
}

sub taxid2division{
  my ($self, $genbank_taxid) = @_;

  if(not $self->taxid2sysname($genbank_taxid)){
    warn "$genbank_taxid is not a valid GenBank taxonomic ID\n";
    return undef;
  }

  my @lineage = $self->lineage($genbank_taxid);
  if($lineage[4] == $metazoa or $lineage[4] == $fungi){
    return $lineage[4];
  }
  if($lineage[3] == $plants){
    return $lineage[3];
  }
  # Look at 3rd level of tree, because 1 and 2 are always root, biota
  elsif($lineage[2] == $bacteria
     or $lineage[2] == $archaea
     or $lineage[2] == $eukarya
     or $lineage[2] == $metagenomes ){
    return $lineage[2];
  }
  elsif($lineage[1] == $viruses){
    return $lineage[1];
  }
  else{
    return 0;
  }
}

sub parent{
  my ($self, $data) = @_;

  if($data =~ /\D/){
    return undef;
  }
  if($data == 0){
    return undef;
  }

  if(not defined $self->{'taxid2parent'}){
      $self->{'taxid2parent'} = {};
      #tie(%{$self->{'taxid2parent'}}, "GDBM_File", "$species_index_path.taxid2parent", O_RDONLY, 0644)
      tie(%{$self->{'taxid2parent'}}, "GDBM_File", "$species_index_path.taxid2parent", O_CREAT|O_WRONLY, 0644)
      #tie(%{$self->{'taxid2parent'}}, "GDBM_File", "$species_index_path.taxid2parent",'rw', 0775)
        or (warn "Cannot open $species_index_path.taxid2parent for reading: $!" and return undef);
  }
  if(!defined ($self->{'taxid2parent'}->{$data} && $login ne "apachebody" && -w "$species_index_path.taxid2parent")){

      print STDERR "$data no parent in local taxonomic db, we will search in NCBI..... " unless $ENV{'MAGPIEQUIET'};

      my $xml = get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$data&report=sgml&mode=xml&email=xdong\@ucalgary.ca";
      print STDERR "$xml\n";
      if($xml =~ /.*?\<ParentTaxId\>(\d+)\<.*/s){
	  sleep 1;  # to not overload the NCBI Web server
          if(!-w "$species_index_path.taxid2parent"){ #can't update the mapping db, just return the result

	      return $1;
          }

	  print STDERR "ncbi taxid=". $self->{'taxid2parent'}->{$1} . ".........";

          if(not defined $self->{'taxid2parent.rw'}){
	      tie(%{$self->{'taxid2parent'}}, "GDBM_File", "$species_index_path.taxid2parent", O_RDWR, 0644)
	      		  or (warn "Cannot open $species_index_path.taxid2parent for read/write: $!" and return undef);
	      $self->{'taxid2parent.rw'} = 1;
          }

          $self->{'taxid2parent'}->{$data} = $1;
	  print STDERR "updated $1\t$data\n";
      }
      print STDERR " done search in NCBI\n" if $ENV{'MAGPIEQUIET'};
  }
  return $self->{'taxid2parent'}->{$data};
}

sub children{
  my ($self, $node) = @_;

  if(not defined $self->{'taxid2children'}){
	$self->{'taxid2children'} = {};
	tie(%{$self->{'taxid2children'}}, "GDBM_File", "$species_index_path.taxid2children", O_RDONLY, 0644)
		or (warn "Cannot open $species_index_path.taxid2children for reading: $!" and return undef);
  }
  return split(/\n/, $_[0]->{'taxid2children'}->{$node});
}

sub all_species{
  my ($self, $node) = @_;
  my @children = ($node);
  my $child;
  for $child ($self->children($node)){
    push @children, $self->all_species($child);
  }
  return grep {$hierarchy{$self->rank($_)} >= $hierarchy{species}} @children;
}

sub rank{
  my ($self, $data) = @_;

  if($data =~ /\D/){
    return undef;
  }

  if(not defined $self->{'taxid2rank'}){
      $self->{'taxid2rank'} = {};
      tie(%{$self->{'taxid2rank'}}, "GDBM_File", "$species_index_path.taxid2rank", O_RDONLY, 0644)
        or (warn "Cannot open $species_index_path.taxid2rank for reading: $!" and return undef);
  }
  return $self->{'taxid2rank'}->{$data};
}

sub gid2taxid{
  my ($self, $gid) = @_;

  if($gid =~ /\D/){
    return undef;
  }

  if(not defined $self->{'gid2taxid'}){
      $self->{'gid2taxid'} = new FileHandle;
      open($self->{'gid2taxid'}, "$species_index_path.gid2taxid")
        or (warn "Cannot open $species_index_path.gid2taxid for reading: $!" and return undef);
  }
  seek $self->{'gid2taxid'}, $gid*4, 0;  #4 byte integers for taxids
  my $data;
  read $self->{'gid2taxid'}, $data, 4;
  my $taxid = unpack("I", $data);
  if($taxid == 0 && $login ne "apachebody" && -w "$species_index_path.gid2taxid"){
    print STDERR "Checking the NCBI for a remapping of GI $gid ..." unless $ENV{'MAGPIEQUIET'};
    my $xml = get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$gid&retmode=xml&email=xdong\@ucalgary.ca" or warn  "Cannot connect to entrez";
    sleep 1;  # to not overload the NCBI Web server
    if($xml =~ /replaced by gi:(\d+)/s){
      print STDERR " -> $1" unless $ENV{'MAGPIEQUIET'};
      $taxid = $self->gid2taxid($1);
    }
    elsif($xml =~ /GBSeq_organism>(.+?)<\//){
      $taxid = $self->name2taxid($1);
    }
    # save the new value to the local db
    if(-w "$species_index_path.gid2taxid"){
      if(not defined $self->{'gid2taxid.rw'}){
        open($self->{'gid2taxid'}, "+<$species_index_path.gid2taxid")
          or (warn "Cannot open $species_index_path.gid2taxid for read/write: $!" and return undef);
#	select(self->{'gid2taxid'});
#	$|=1;
       $self->{'gid2taxid.rw'} = 1;
      }
      seek $self->{'gid2taxid'}, $gid*4, 0;
      print {$self->{'gid2taxid'}} pack("I", $taxid);
      print STDERR " taxid: $taxid\n" unless $ENV{'MAGPIEQUIET'};

    }
  }
  return $taxid;
}

sub taxid2sysname{
  my ($self, $data) = @_;

  if($data =~ /\D/){
    return undef;
  }

  if(not defined $self->{'taxid2sysname'}){
      $self->{'taxid2sysname'} = {};
      #tie(%{$self->{'taxid2sysname'}}, "GDBM_File", "$species_index_path.taxid2sysname", O_RDONLY, 0644)
      tie(%{$self->{'taxid2sysname'}}, "GDBM_File", "$species_index_path.taxid2sysname", O_CREAT|O_WRONLY, 0644)
        or (warn "Cannot open $species_index_path.taxid2sysname for reading: $!" and return undef);
  }
  my $sysname = $self->{'taxid2sysname'}->{$data};

  if($sysname =~ /ParentTaxId/ || length($sysname) == 0){

      my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$data&report=sgml&mode=xml&email=xdong\@ucalgary.ca";
      my $xml = get $url;

      print STDERR "$url\n\n$xml";
      if($xml =~ /\<ScientificName\>?(\S.*?)\<\/ScientificName\>/s){
	  if(length($1) > 0){
	      print STDERR "xiaoli data=$data********************$1***********end xiaoli\n";
	      $self->{'taxid2sysname'}->{$data} = $1;
	      $sysname = $1;
	      return $sysname;
	  }
      }
  }
  return $sysname;
}

# If name is ambiguous, first one encountered is returned.
# To disambiguate, pass in the parent's name or taxid too.
sub name2taxid{
    my ($self, $data, $parent) = @_;
    #print STDERR "$data, $parent\n";
    if(not defined $self->{'name2taxid'}){
	$self->{'name2taxid'} = {};
	#tie(%{$self->{'name2taxid'}}, "GDBM_File", "$species_index_path.name2taxid", O_RDONLY, 0644)
        #or (warn "Cannot open $name2taxid for reading: $!" and return undef);
	tie(%{$self->{'name2taxid'}}, "GDBM_File", "$species_index_path.name2taxid", O_CREAT|O_WRONLY, 0644)
	    or (warn "Cannot open $name2taxid for reading: $!" and return undef);


    }

    my @taxids = split /\n/, $self->{'name2taxid'}->{lc($data)};
    if(@taxids < 1){
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=$data&mode=xml&email=xdong\@ucalgary.ca&report=sgml";
	my $xml = get $url;
	#print STDERR "$url\n\n$xml";
	if($xml =~ /\<IdList\>(.*)\<\/IdList\>/s){

	    my @idlist = split(/\n/, $1);
	    foreach my $idline (@idlist){
		if(my ($id) = $idline =~ /\<Id\>(\d+)\<\/Id\>/s){
		    print STDERR "idline=$id\n";
		    if(defined $self->{'name2taxid'}->{lc($data)}){
			$self->{'name2taxid'}->{lc($data)} .= "\n$id";
		    }
		    else{
			$self->{'name2taxid'}->{lc($data)} = $id;
		    }
		}

	    }
	}
	sleep 1;  # to not overload the NCBI Web server
    }


    #print STDERR "magpie $data********** @taxids\n";
    #if(@taxids == 1 or not $parent or $parent eq $data){
	#print STDERR "$data********** $taxids[0]\n";
#	return $taxids[0];
 #   }
    if(@taxids == 1){
	#print STDERR "$data********** $taxids[0]\n";

	return $taxids[0];
    }

    else{
	# The parent could be specified as text or a taxid
	my ($parent_taxid) = $parent =~ /^\d+$/ ? $parent : $self->{'name2taxid'}->{lc($parent)};

#print STDERR "parent=$parent\t$parent_taxid\n";

	my $shortest_lineage;
	my $shortest_lineage_taxid;
	my $shortest_distance = 999;

	for my $taxid (@taxids){

	    my $lineage_distance = 0;

	    for my $lineage_taxid (reverse($self->lineage($taxid))){
		$lineage_distance++;
		if($parent_taxid == $lineage_taxid){
		    return $lineage_taxid if $lineage_distance == 1; # direct parent
		    last;
		}
	    }
	    if(not defined $shortest_distance or $lineage_distance < $shortest_distance){
		#print STDERR "&&&&&&&&&&&&&&& $taxid=$lineage_distance;\n";
		$shortest_distance = $lineage_distance;
		$shortest_lineage_taxid = $taxid;
	    }
	}
	# no parent match, return first one anyways
	return $shortest_lineage_taxid;
    }
}




sub TaxFileDir{
  return $ENV{'MAGPIEHOME'}."/Tables/Taxonomy";
}

sub updateTaxonomyIndices{
  my $file = TaxFileDir()."/names.dmp";

  open(TAXDB, $file)
    or (warn "Cannot open $file for reading: $!\n" and return undef);
  my (%taxid2sysname, %name2taxid);

  # Delete the old files
  unlink "$species_index_path.taxid2sysname"
    or (warn "Cannot delete old $species_index_path.taxid2sysname: $!" and return undef);
  unlink "$species_index_path.name2taxid"
    or (warn "Cannot delete old $species_index_path.name2taxid: $!" and return undef);

  # English names mapped to taxonomic ids
  tie(%taxid2sysname, "GDBM_File", "$species_index_path.taxid2sysname", O_RDWR|O_CREAT, 0644)
    or (warn "Cannot open $species_index_path.taxid2sysname for writing: $!" and return undef);
  tie(%name2taxid, "GDBM_File", "$species_index_path.name2taxid", O_RDWR|O_CREAT, 0644)
    or (warn "Cannot open $species_index_path.name2taxid for writing: $!" and return undef);
  while($_ = <TAXDB>){
    my @fields = split /\t\|\t/;
    if($fields[3] eq "scientific name\t|\n"){
      $taxid2sysname{$fields[0]} = $fields[1];
      #print STDERR lc($fields[1]) , " (canonical)\t", $., "\n";
      if(defined $name2taxid{lc($fields[1])}){
          $name2taxid{lc($fields[1])} .= "\n$fields[0]";
      }
      else{
          $name2taxid{lc($fields[1])} = $fields[0];
      }

      # There are cases where a name might be assigned to more than one classfication
      # level: in this case also add a bare name entry for the one with (class) (phylum), etc. at the end
      if($fields[1] =~ /\(class\)/i){
	  my $bare_name = lc($fields[1]);
	  $bare_name =~ s/\s+\(class\)\s*$//; # strip it
	  #print STDERR "Adding alias $bare_name for $fields[1]\n";
	  if(defined $name2taxid{$bare_name}){
	      $name2taxid{$bare_name} .= "\n$fields[0]";
	  }
	  else{
	      $name2taxid{$bare_name} = "$fields[0]";
	  }
      }
    }
    elsif($#fields == 3){
      my $name = lc($fields[1]);
      $name =~ s/^(["'])(.+?)\1.*$/$2/; #Crap after names in quotes
      $name =~ s/^(.*) \S+ (?:and \S+|et al\.) (?:1[6-9]\d\d|200\d).*/$1/g;
      $name =~ s/\s+\((?:ex |Approved ).*?\)//g;
      $name =~ s/[^A-Za-z0-9\)\.]$//;
      $name =~ s/\((?:\s*|sic)\)//g;
      $name =~ s/\([^\)]+$//;
      $name =~ s/\s+/ /g;
      #print STDERR $name , "\t", $., "\n";
      if(defined $name2taxid{$name} and $name2taxid{$name} !~ /\b$fields[0]\b/){
          $name2taxid{$name} .= "\n$fields[0]";
      }
      else{
          $name2taxid{$name} = $fields[0];
      }
    }
    else{
      print STDERR "Line $. of $file unparseable";
    }
  }
  untie %taxid2sysname;
  untie %name2taxid;

  # Hierarchy info for names
  $file = TaxFileDir()."/nodes.dmp";
  open(TAXDB, $file)
    or (warn "Cannot open $file for reading: $!\n" and return undef);
  my (%taxid2parent, %taxid2children, %taxid2rank);

  # Delete the old files
  if(-e "$species_index_path.taxid2parent"){
    unlink "$species_index_path.taxid2parent"
      or (warn "Cannot delete old $species_index_path.taxid2parent: $!" and return undef);
  }
  if(-e "$species_index_path.taxid2children"){
    unlink "$species_index_path.taxid2children"
      or (warn "Cannot delete old $species_index_path.taxid2children: $!" and return undef);
  }
  if(-e "$species_index_path.taxid2rank"){
    unlink "$species_index_path.taxid2rank"
      or (warn "Cannot delete old $species_index_path.taxid2rank: $!" and return undef);
  }

  tie(%taxid2parent, "GDBM_File", "$species_index_path.taxid2parent", O_RDWR|O_CREAT, 0644)
    or (warn "Cannot open $species_index_path.taxid2parent for writing: $!" and return undef);
  tie(%taxid2children, "GDBM_File", "$species_index_path.taxid2children", O_RDWR|O_CREAT, 0644)
    or (warn "Cannot open $species_index_path.taxid2children for writing: $!" and return undef);
  tie(%taxid2rank, "GDBM_File", "$species_index_path.taxid2rank", O_RDWR|O_CREAT, 0644)
    or (warn "Cannot open $species_index_path.taxid2rank for writing: $!" and return undef);
  while($_ = <TAXDB>){
    if(/^(\d+)\t\|\t(\d+)\t\|\t(\S+)/){
      $taxid2parent{$1} = $2;
      $taxid2rank{$1} = $3;
      if(not defined $taxid2children{$2}){
        $taxid2children{$2} = $1
      }
      else{
        $taxid2children{$2} .= "\n$1";
      }
    }
    else{
      print STDERR "Line $. of $file unparseable";
    }
  }
  untie %taxid2parent;
  untie %taxid2children;
  untie %taxid2rank;

  # Map gi numbers to species
  my $gid2taxid = new FileHandle;
  open($gid2taxid, ">$species_index_path.gid2taxid")
     or (warn "Cannot open $species_index_path.gid2taxid for writing: $!" and return undef);
  for $file (qw(gi_taxid_nucl.dmp gi_taxid_prot.dmp)){
    my $datafile = TaxFileDir()."/$file";
    open(TAXDB, $datafile)
      or (warn "Cannot open $datafile for reading: $!\n" and return undef);
    while($_ = <TAXDB>){
      chomp;
      my($g, $t) = split /\t/;
      next if $t eq "0"; #Unmapped
      seek $gid2taxid, $g*4, 0;  #4 byte integers for taxids
      print $gid2taxid pack("I", $t);
    }
  }
  close(TAXDB);
  close($gid2taxid);

#TODO: chmod to 0644 of the database

  return 1;
}

sub isSpecies{
  my ($self, $data);
  if($#_){  #Was called as object instance
    ($self, $data) = @_;
  }
  else{
    $self = {};
    $data = $_[0];
  }

  if(length $data < 3 or $data eq "beta"){ #Crap and common word species shortcut
    return 0;
  }

  #Obvious abbreviation cases, H. sapiens or Anas sp., assuming you're only feeding description lines!
  if($data =~ /^[A-Z]\.\s+[a-z]{3,}$/ or $data =~ /^[A-Z][a-z]+\s+sp\.$/){
    return 1;
  }

  #See if it's a systematic or common/alternative name
  if(not defined $self->{'name2taxid'}){
    $self->{'name2taxid'} = {};
    tie(%{$self->{'name2taxid'}}, "GDBM_File", "$species_index_path.name2taxid", O_RDONLY, 0644)
      or die "Cannot open $species_index_path.name2taxid for reading: $!";
  }
  if(defined $self->{'name2taxid'}->{lc($data)}){
    return 1;
  }
  else{
    return 0;
  }
}

1;
