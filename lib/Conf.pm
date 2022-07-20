#This library contains subroutines for reading in the configuration files of 
#a project.
#Paul Gordon
#
# 	$Id: Conf.pm,v 1.17 2010/11/05 16:18:40 gordonp Exp $	
package Conf;

#require Exporter;
#@EXPORT = qw(project);
require 'messages.pl';
use strict qw(vars);
use Constants;
use Tool;

sub Conf::magpiehome(){
  my $magpiehome = $ENV{'MAGPIEHOME'} 
    || main'fatal("The environment variable MAGPIEHOME is not defined");

  if($magpiehome =~ /^([^;\|\`\n]+)$/){
    $magpiehome = $1;
  }

  main'fatal("The MAGPIEHOME '$magpiehome' does not exist") 
    unless -e $magpiehome;

  main'fatal("The MAGPIEHOME '$magpiehome' is not a directory or link") 
    unless -d $magpiehome || -l $magpiehome;
  
#  main'fatal("The MAGPIEHOME '$magpiehome' doesn't appear to contain a full",
#	"MAGPIE distribution") unless -e "$magpiehome/rc" && 
#	                              -e "$magpiehome/bin" && 
#	                              -e "$magpiehome/Config" && 
#	                              -e "$magpiehome/lib";

  return $magpiehome;
}

sub www_config{
  my %conf;
  my $file = &magpiehome ."/Config/WWWCONFIG";

  if(not open(WWWCFG, $file)){
    warn "Cannot open $file for reading: $!\n";
    return undef;
  }
  while(<WWWCFG>){
    next if /^\s*(?:\#.*)?$/;
    chomp;
    if(/^(\S+)\s*:\s*(\S+)/){
      $conf{$1} = $2; 
    }
    else{
      warn "Malformatted line $. in $file: should be foo:bar, found $_\n";
    }
  }
  close(WWWCFG);

  return \%conf;
}

sub public_groups{
  my ($project) = @_;
  my %conf;
  my $file = $project->config->{'path'}."/Config/PUBLICCONFIG";
  
  if(not open(PUBCFG, $file)){
    warn "Cannot open $file for reading: $!\n";
    return undef;
  }
  while(<PUBCFG>){
    next if /^\s*(?:\#.*)?$/;
    if(/^(\S+)\s+(\S+)\s*$/){
      $conf{$1} = $2;
    }
    else{
      warn "Unrecognized line (#$.) in $file: ignoring\n"; 
    }
  }
  close(PUBCFG);
  return \%conf;
}

sub web_align{
  my ($scope) = @_; 
  my (%conf, $file);

  if($scope == $Constants::private){
    $file = magpiehome()."/cgi-bin/private/.alignrc";
  }
  else{
    $file = magpiehome()."/cgi-bin/public/.alignrc";
  }
  if(not open(ALIGNRC, $file)){
    warn "Cannot open $file for reading: $!\n";
    return undef;
  }

  while(<ALIGNRC>){
    next if /^\s*(?:\#.*)?$/;
    my @fields = split /\t/, $_;
    my $field;
    foreach $field (@fields){
      if($field eq "AA"){
         $field = $Constants::protein;
      }
      elsif($field eq "DNA"){
         $field = $Constants::dna;
      }
      elsif($field eq "DNAs"){
         $field = $Constants::dnascope;
      }
      elsif($field eq "ORFs"){
         $field = $Constants::orfscope;
      }
      elsif($field eq "AAs"){
         $field = $Constants::aascope;
      }
    } 
    if($#fields == 5){
      if(not exists $conf{$fields[0]}){
        $conf{$fields[0]} = {}; 
      }
      if(not exists $conf{$fields[0]}->{$fields[1]}){
	$conf{$fields[0]}->{$fields[1]} = {};
      }	
      if(not exists $conf{$fields[0]}->{$fields[1]}->{$fields[2]}){
        $conf{$fields[0]}->{$fields[1]}->{$fields[2]} = {};
      } 
      if(not exists $conf{$fields[0]}->{$fields[1]}->{$fields[2]}->{$fields[3]}){
        $conf{$fields[0]}->{$fields[1]}->{$fields[2]}->{$fields[3]} = [@fields[4..5]];
      } 
    }
  }

  return \%conf; 
}

sub pbs_config{
  my ($project) = @_;

  my %conf;
  my $file = $project->config->{'path'}."/Config/PBSCONFIG";
  if(not -e $file){
    $file = $ENV{'MAGPIEHOME'}."/Config/PBSCONFIG";
  }

  if(not -e $file){
    return undef;
  }

  open(CFG, $file)
    or die "Cannot open $file for reading: $!\n";

  while(<CFG>){
    if(/^(access|time|qsub|subargs)=(.*)$/){
      $conf{$1} = $2;
    }
  }

  return \%conf;  
}

#Load all of the special words that the BioWords module might need
sub keyWords{
  my ($project, $biowords) = @_;

  my $dir;
  for $dir ($project->config->{'path'}, $ENV{'MAGPIEHOME'}){
    if(-e "$dir/Tables/Words/stops"){
      open(STOPFILE, "$dir/Tables/Words/stops")
        or warn "Cannot open $dir/Tables/Words/stops for reading: $!";
      $biowords->addStopWordsFromFile(\*STOPFILE);
      close(STOPFILE);
    }
    if(-e "$dir/Tables/Words/kills"){
      open(KILLFILE, "$dir/Tables/Words/kills")
        or warn "Cannot open $dir/Tables/Words/kills for reading: $!";
      $biowords->addKillPatternsFromFile(\*KILLFILE);
      close(KILLFILE);
    }
  }
}

# Returns a hash table with uppercase codon as key, array with two elements,
# codon frequency proportions vs. all codons, and vs. codons for this amino acid
sub codon_usage{
  my ($project) = @_;
  my %conf;
  my $file = magpiehome()."/Tables/GeneSkipper/".$project->config->{'codontable'
}.".cod";

  if(not open(CUT, $file)){
     warn "Cannot open $file for reading: $!\n";
     return undef;
  }

  while(<CUT>){
    next if /^(#|\s*$)/;
    if(/^\S+\s+([AGCT]{3})\s+\S+\s+([0-9\.]+)\s+([0-9\.]+)/i){
      $conf{uc($1)} = [$2/1000, $3];
    }
  }
  close(CUT);

  return \%conf;
}

sub translation_table{
  my ($project) = @_; 
  my %conf;
  my $file = magpiehome()."/Tables/GeneSkipper/".$project->config->{'codontable'}.".cod";

  my $oldsep = $/;
  $/ = "\n";
  if(-e $file)
  {
      open(CUT, $file)
        || warn "Cannot open codon usage table $file for reading: $!\n";
      while(<CUT>)
         {if(/^(\S+)\s+([ACGT]{3})/i){$conf{uc $2} = $1;}}
      close(CUT);
  }
  else{
    %conf = (
'TTT', 'Phe', 'TCT', 'Ser', 'TAT', 'Tyr', 'TGT', 'Cys',
'TTC', 'Phe', 'TCC', 'Ser', 'TAC', 'Tyr', 'TGC', 'Cys',
'TTA', 'Leu', 'TCA', 'Ser', 'TAA', 'Stop','TGA', 'Stop',
'TTG', 'Leu', 'TCG', 'Ser', 'TAG', 'Stop','TGG', 'Trp',
'CTT', 'Leu', 'CCT', 'Pro', 'CAT', 'His', 'CGT', 'Arg',
'CTC', 'Leu', 'CCC', 'Pro', 'CAC', 'His', 'CGC', 'Arg',
'CTA', 'Leu', 'CCA', 'Pro', 'CAA', 'Gln', 'CGA', 'Arg',
'CTG', 'Leu', 'CCG', 'Pro', 'CAG', 'Gln', 'CGG', 'Arg',
'ATT', 'Ile', 'ACT', 'Thr', 'AAT', 'Asn', 'AGT', 'Ser',
'ATC', 'Ile', 'ACC', 'Thr', 'AAC', 'Asn', 'AGC', 'Ser',
'ATA', 'Ile', 'ACA', 'Thr', 'AAA', 'Lys', 'AGA', 'Arg',
'ATG', 'Met', 'ACG', 'Thr', 'AAG', 'Lys', 'AGG', 'Arg',
'GTT', 'Val', 'GCT', 'Ala', 'GAT', 'Asp', 'GGT', 'Gly',
'GTC', 'Val', 'GCC', 'Ala', 'GAC', 'Asp', 'GGC', 'Gly',
'GTA', 'Val', 'GCA', 'Ala', 'GAA', 'Glu', 'GGA', 'Gly',
'GTG', 'Val', 'GCG', 'Ala', 'GAG', 'Glu', 'GGG', 'Gly');
  }

  $/ = $oldsep;

  return \%conf;
}

sub coverage_config{
  my ($project) = @_;

  my %conf;
  my $file = $project->config->{'path'}."/Config/COVERAGECONFIG";
  if(not -e $file){
     $file = $ENV{MAGPIEHOME}."/Config/COVERAGECONFIG";
  }
  open(CFG, $file) or die "Cannot open $file for reading: $!\n";

  while(<CFG>){
    next if /^\s*(#|$)/;
    chomp;
    my @fields = split /\s+/;
    if($#fields == 3){  
      $conf{$fields[0]} = [@fields[1..3]]; 
    }
    elsif(/^(LEVEL[123])\s*=\s*([0-9\.]+)\s*$/){
      $conf{$1} = $2;
    }
    elsif(/^TYPE:/){
      # Ignore old format
    }
    else{
      warn "Malformed line ($file #$.): expected four spaced fields\n";
    }
  }
  close(CFG);

  $conf{LEVEL1} ||= 1;

  return \%conf;
}

sub score_config{
  my ($project) = @_;

  my %conf;
  my $file = $project->config->{'path'}."/Config/SCORECONFIG";
  if(not -e $file){
     $file = $ENV{MAGPIEHOME}."/Config/SCORECONFIG";
  }
  open(CFG, $file) or die "Cannot open $file for reading: $!\n";

  my $oldsep = $/;
  $/ = "\n";
  while(<CFG>){
    next if /^\s*(?:#|$)/;
    chomp;
    if(/^(\S+)\s+combine\s+by\s+and/){
      $conf{$1}->{and} = 1;
    }
    else{
      if(/(\S+)\s+([1234])\s+(\S+)\s+(\S+)\s+(".*?"|\S+)/){
        if(not defined $conf{$1}){
          $conf{$1} = {};
        }
        if(not defined $conf{$1}->{$2}){
           $conf{$1}->{$2} = "\$facts->{'$3'} $4 $5";
        }
        elsif($conf{$1}->{and}){
           $conf{$1}->{$2} .= " and \$facts->{'$3'} $4 $5";
        }
        else{
           $conf{$1}->{$2} .= " or \$facts->{'$3'} $4 $5";
        }
      }
      elsif(/(\S+)\s+([1234])\s+all/){
        $conf{$1}->{$2} = 1;
      }
      else{
        warn "Malformatted line ($file #$.): expected 'tool level measure op value'"; 
      }
    }
  }
  $/ = $oldsep;

  return \%conf;
}

sub set_project($$){
  my($projectname, $settings_ref) = @_;
  my $file = magpiehome()."/rc/.magpierc.$projectname";

  if($file =~ /^([A-Za-z0-9\/_\-\.]+)$/){ 
      $file = $1;
  }
  else{ die "Filename: $file contains illegal characters.";}  
  if(not open(CFG, ">$file")){
    warn "Cannot open $file for writing: $!\n";
    return 0;
  }
  my $key;
  for $key (sort keys %$settings_ref){
    print CFG "$projectname $key:", $settings_ref->{$key}, "\n";
  }
  close(CFG);

  chmod 0644, $file;

  return 1;
}

sub Conf::project($$){
  my($project, $settings_ref) = @_;
  my $magpiehome = magpiehome();
  my %new_setting;
  my %settings =   ('path' => undef, 
		    'mailqueue' => undef,
		    'input' => undef,
                    'gbrowse' => '',
		    'notify' => "root\@localhost",
		    'mailaddr' => '',
		    'codontable' => 'ecoli',
		    'startcodons' => '"ATG","GTG","TTG"',
		    'stopcodons' => '"TAG","TGA","TAA"',
		    'organism' => 'unknown',
		    'overlapping_orf_cutoff' => '1.0',
		    'promoter_search' => 'no',
		    'shinedel_search' => 'no',
		    'href' => '',
		    'cgi-bin' => undef,
		    'hrefbase' => undef,
		    'httpbasedir' => $magpiehome."/project_links/$project",
		    'rare_max' => 25,
		    'ag_min' => 50,
		    'cu_min' => 37,
		    'gc_diff_min' => 10,
		    'db' => '',
		    'mach' => '',
		    'min_intron' => '',
		    'whole_genome' => 'no',
		    'augustus_model' => '',
		    'glimmermodel' => '',
                    'glimmerlength' => '',
                    'min_interorf' => 30,
		    'ec_lookup' => 0,
		    'ec_lookup_local' => 'no',
		    'vector' => [],
                    'orfid' => '',
                    'dynamic' => 0,
                    'queuefile' => '',
                    'bae_division' => '',
                    'annotation_equiv_pct_query' => 95,
                    'annotation_equiv_pct_subject' => 95,
                    'annotation_equiv_pct_sim' => 95,
                    'board' => 0,
                    'pbs' => 0,
                    'indexed' => 0,
                    'taxid' => 0);

  #delete any previous settings
  for my $key (keys %$settings_ref){
    delete $$settings_ref{$key};
  }
  #set defaults
  for my $key (keys %settings){
    $$settings_ref{$key} = $settings{$key};
  }

  my $rcfile = "$magpiehome/rc/.magpierc.$project";
  main'fatal("The project '$project' does not exist.") unless -e $rcfile;

  open (CONFIG, $rcfile) 
    or main'fatal("Cannot open project configuration file '$rcfile' for reading: $!");
  my $oldsep = $/;
  $/ = "\n";
  while(<CONFIG>){
    next if /^\s*(?:\#.*)?$/;
#    if(my ($setting, $value) = /^\s*$project\s+(\S+?)\s*:\s*(\S.*?)\s*$/){
#      main'fatal("Found unrecognized setting '$setting' in $rcfile line $.")
#	unless exists $settings{$setting}; 
#      if($setting eq 'vector'){
#	push @{$$settings_ref{$setting}}, $value;
#      }
#      elsif(defined $new_setting{$setting}){
#	main'warning("Duplicate value in $rcfile for '$setting' (", $value,
#		"): using first definition (", $$settings_ref{$setting},")");
#      }
#      else{
#	$new_setting{$setting} = $$settings_ref{$setting} = $value;
#      }
#    }
#    else{
#      main'warning("Unparseable line ('$rcfile' #$.), expected ",
#	      "\"$project setting=value\"");
#    }
    if(/^\s*$project\s+(\S+?)\s*:\s*(\S.*?)\s*$/){
      main'fatal("Found unrecognized setting '$1' in $rcfile line $.")
	unless exists $settings{$1}; 
      if($1 eq 'vector'){
	push @{$$settings_ref{$1}}, $2;
      }
      elsif(defined $new_setting{$1}){
	main'warning("Duplicate value in $rcfile for '$1' (", $2,
		"): using first definition (", $$settings_ref{$1},")");
      }
      else{
	$new_setting{$1} = $$settings_ref{$1} = $2;
      }
    }
    else{
      main'warning("Unparseable line ('$rcfile' #$.), expected ",
	      "\"$project setting:value\"");
    }
  }
  close (CONFIG);
  $/ = $oldsep;

  my @missing;
  for my $key (keys %$settings_ref){
    push @missing, $key if not defined $$settings_ref{$key};
    if($$settings_ref{$key} =~ /^yes$/i){$$settings_ref{$key} = 1}
    elsif($$settings_ref{$key} =~ /^no$/i){$$settings_ref{$key} = 0}
  }
  main'fatal("Required setting",($#missing?'s':'')," missing in '$rcfile':", 
	join(",", @missing)) if @missing;
}

sub graphics(){
  my ($project, $settings_ref, $project_settings_ref) = @_;
  my $context = "";
  my %script_setting = ();

  my ($script) = $0 =~ /([^\/]+)$/;
  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &read_project_config($project, $project_settings_ref) 
  }

  my $file = $$project_settings_ref{'path'}."/Config/GRAPHCONFIG";
  # Use the installation config if no project config exists
  if(!-e $file){  
    $file = &magpiehome()."/Config/GRAPHCONFIG";
  }
  open(CONFIGFH, $file) 
    or main'fatal("Cannot open graphics config file '$file' for reading: $!");

  my $line;
  while($line = <CONFIGFH>){
    next if $line =~ /^\s*(?:\#.*)?$/;
    if($line =~ /^start\s+($script|global)/i){
      $context = $1;
      while(<CONFIGFH>){
        next if /^\s*(?:\#.*)?$/;
        if(/^\s*end\s+\S+/i){
          $context = '';
          last;
        }
        if(my($setting, $value) = /^\s*([^\#]\S*?)\s*=\s*(.*?)\s*$/){
          if(!exists $script_setting{$setting}){
	    if($context eq "global" and exists $$settings_ref{$setting}){
	      main'warning("Global setting '$setting' redefined ('$file' #$.), ",
		      "using first value '", $$settings_ref{$setting}, "'")
	    }
	    else{
	      $$settings_ref{$setting} = $value;
	    }
	  }
          if($context eq $script){
            main'warning("Setting '$setting' for $script redefined ('$file' ",
		    "#$.), using first value '", $$settings_ref{$setting}, "'")
	      if $script_setting{$setting}++;
	  }
        }
        else{
          main'warning("Improperly formatted line ('$file' #$.), expected 'key = value'");
        }
      }
    }
  }
  close(CONFIGFH);
}

sub read_color_config($$){
  my($project, $settings_ref, $project_settings_ref) = @_; 

  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &read_project_config($project, $project_settings_ref) 
  }

  my $file = $$project_settings_ref{'path'}."/Config/COLORCONFIG";
  # Use the installation config if no project config exists
  if(!-e $file){  
    $file = &magpiehome()."/Config/COLORCONFIG";
  }
  open(CONFIGFH, $file)
    or main'fatal("Cannot open color config file '$file' for reading: $!");

  while(<CONFIGFH>){
    next if /^\s*(?:\#.*)?$/;
    if(my ($tool, $r, $g, $b) = /^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/){
      if(exists $$settings_ref{$tool}){
	main'warning("Tool color redefine ('$file' #$.) for $tool, using first ",
	        "definition R:", $$settings_ref{$tool}[0], 
		          " G:", $$settings_ref{$tool}[1], 
		          " B:", $$settings_ref{$tool}[2]);
      }
      else{
	$$settings_ref{$tool} = [$r, $g, $b];
      }
    }
  }
  close(CONFIGFH);
}

sub read_priority_config{
  my($project, $project_settings_ref) = @_;
  my (%run_priority, %run_hosts, %run_steps);
  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &read_project_config($project, $project_settings_ref)
  }
  
  my $priorityfile = $$project_settings_ref{'path'}."/Config/PRIORITYCONFIG";
  # Use the installation config if no project config exists
  if(!-e $priorityfile){
    $priorityfile = &magpiehome()."/Config/PRIORITYCONFIG";
  } 
  open(CONFIGFH, $priorityfile)
    or main'fatal("Cannot open tool priority file '$priorityfile' for reading: $!");
  while(<CONFIGFH>){
    next if /^\s*(?:\#|$)/;
    if(/^\s*(\S+)\s+(\d+)\s*(\S*)\s*(\d*)\s*$/){
      $run_priority{$1} = $2;
      $run_hosts{$1} = $3;
      $run_steps{$1} = $4;
    }
  }
  close(CONFIGFH);
  return(\%run_priority, \%run_hosts, \%run_steps);
}

sub read_tool_config{
  my($project, $tool_settings_ref, $project_settings_ref) = @_;
  my($aa, $dna, $hmm, $msa) = qw(AA DNA HMM MSA);
  my($text, $fact) = qw(text fact);
  my %simmap = qw(prosim protein
                  dnasim dna
                  motif  motif
                  model  model
                  trna   trna
                  repeat repeat
                  orf orf);

  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &read_project_config($project, $project_settings_ref) 
  }

  my ($run_priority, $run_hosts, $run_steps) = &read_priority_config($project, $project_settings_ref);
  my $part;
  for $part (qw(TOOLCONFIG InterORFCONFIG)){
  
  my $file = $$project_settings_ref{'path'}."/Config/$part";
  # Use the installation config if no project config exists
  if(!-e $file){  
    $file = &magpiehome()."/Config/$part";
  }  

  open(CONFIGFH, $file)
    or main'fatal("Cannot open tool config file '$file' for reading: $!");

  my $priority = 0;
  $/ = "\n";
  while(<CONFIGFH>){
    $priority++;
    next if /^\s*(?:\#|$)/;
    my $scope = ($part eq "TOOLCONFIG" ? $Constants::orfscope : $Constants::interorfscope);
    if(my ($name, $in_mol, $runtype, $datatype, $tooltype, $parser, $texttype) = 
       /^\s*(\S+)\s+         # tool name
       ($aa|$dna|$hmm|$msa)\s+         # input type for program query
       (\S+)\s+ # local, e-mail, or http server run tool
       ($fact|$text)\s+        # digested facts are in the feature files
       (\S+)\s+              # tool type (e.g. dnasim)
       (\S+)\s*(toxml|xml)?\s*$/xio){           # tool output format (has corresponding 
                             #  digest_format script in magpie/bin)



          if($in_mol eq $hmm){
            #print "$project, $name, $in_mol, $runtype\n";
          }
      if($datatype eq "fact"){
        if(exists $tool_settings_ref->{$name}){
          if($tool_settings_ref->{$name}->scope == $scope){
            main'warning("Redefinition for fact tool '$name' in tool config file ",
                  "('$file' #$.), will use first definition");
          }
          else{
            $scope = $Constants::allscope;
          }
        }
        else{
          if(not $simmap{$tooltype}){
            die "At line $. of $file: $tooltype is not a valid type, must be one of ",
                join(", ", keys %simmap);
          }


          my $tool = new Tool($project, $name, ($in_mol eq $aa?$Constants::protein : $in_mol eq $hmm? $Constants::hmm : $in_mol eq $msa ? $Constants::msa : $Constants::dna),
                                               $runtype,
                                               ($datatype eq $fact?$Constants::facttool : $Constants::texttool),
                                               ${"Constants::".$simmap{$tooltype}},
                                               $scope,
                                               $priority,
                                               $run_priority->{$name},
                                               $run_hosts->{$name},
                                               $run_steps->{$name},
                                               $parser, $texttype);

#print STDERR $tool->run_hosts. " 00000000000" . $tool->parser . " 000000000" . $tool->name . "\n";
          if(not $tool){
            die "Error creating Tool object corresponding to line $. of $file";
          }
          $tool_settings_ref->{$name} = $tool;

        }
      }
    }
    elsif(!/^\s*\S+\s+[a-z]$/){
      main'warning("Unrecognized line in tool config file ('$file' #$.), expected ",
              "'toolname (AA|DNA|HMM) (local|e-mail\@addr) (text|fact) ",
              "tool_type tool_format'\n or 'tool_type one_letter_abbrev'");
    }
  }


  close(CONFIGFH);
  } #end foreach $part

}

sub read_prefix_config{
  my($project, $file_settings_ref, $project_settings_ref) = @_;

  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &read_project_config($project, $project_settings_ref)
  }

  my $key;
  for $key (keys %$file_settings_ref){
    delete $$file_settings_ref{$key}
  }

  my $cfgfile = $$project_settings_ref{"path"}."/Config/FILECONFIG";
  
  open(CFG, $cfgfile)
    or main'fatal("Cannot open $cfgfile for reading: $!");
  while(<CFG>){
    next if /^\s*(?:\#|$)/;
    if(/\s*(\S+)\s+(\S+)\s*$/){
       if(exists $$file_settings_ref{$1}){
          main'warning("Prefix for $1 redefined ($cfgfile #$.), ignoring");
       }
       else{
	  $$file_settings_ref{$1} = $2;
       }
    }
    else{
       main'warning("Unrecognized line in $cfgfile (#$.), expecting 'group prefix'");
    }
  }
  close(CFG);

}

sub add_prefix{
  my ($project, $groupname, $prefix) = @_;

  my $cfile = $project->config->{'path'}."/Config/FILECONFIG";
  open(CFILE, $cfile) or return undef;
  my @lines = <CFILE>;

# used to pass taint checking
  if($cfile =~ /^([A-Za-z0-9\/_\-\.]+)$/){ 
      $cfile = $1;
  }
  else{
      die "The file: $cfile contains illegal charaters.";
  }

  open(CFILE, ">$cfile") or return undef;
  #Remove any previous file association
  @lines = grep {!/^$groupname\s/} @lines;
  print CFILE @lines, "$groupname\t$prefix\n";
  close(CFILE);
  return 1;
}

sub read_offset_config{
  my($project, $offsets_ref, $project_settings_ref) = @_;

  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &read_project_config($project, $project_settings_ref)
  }

  my $cfgfile = $$project_settings_ref{"path"}."/Config/OFFSETCONFIG";
  open(CFG, $cfgfile)
    or main'fatal("Cannot open $cfgfile for reading: $!");
  while(<CFG>){
    next if /^\s*(?:\#|$)/;
    if(/^\s*(\S+)\s+(\S+)\s+(-?\d+)\s*$/){
       if(not defined $$offsets_ref{$1}){
          $$offsets_ref{$1} = {};
       }
       if(exists $$offsets_ref{$1}->{$2}){
          main'warning("Offset for $1 sequence $2 redefined ($cfgfile #$.), ignoring");
       }
       else{
          $$offsets_ref{$1}->{$2} = $3;
       }
    }
    else{
       main'warning("Unrecognized line in $cfgfile (#$.), expecting 'group sequence offset'");
    }
  }
  close(CFG);

}

sub add_offset{
  my($project, $groupname, $offset) = @_;

  my $cfgfile = $project->config->{"path"}."/Config/OFFSETCONFIG";
  open(CFG, $cfgfile)
    or return undef;
  my @lines = <CFG>;
  #Remove any previous offset for this group
  @lines = grep {!/^$groupname\s/} @lines;
  open(CFG, ">$cfgfile")
    or return undef;
  print @lines, "$groupname\t$offset\n";
  close(CFG);
  return 1;
}

sub read_state_config{
  my($project, $state_settings_ref, $project_settings_ref) = @_;
  
  my $tmpln = $/;
  $/ = "\n";

  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &read_project_config($project, $project_settings_ref) 
  }
  
  my $cfgfile = $$project_settings_ref{"path"}."/Config/STATECONFIG";
  # Use the installation config if no project config exists
  if(!-e $cfgfile){  
    $cfgfile = &magpiehome()."/Config/STATECONFIG";
  }
  open(CFG, $cfgfile)
    or main'fatal("Cannot open $cfgfile for reading: $!");
  my $statename = "";
  while(<CFG>){
    next if /^\s*(?:\#|$)/;
    if(/^(\S+)\s*$/){

      if($1 eq "proteinAA"){
        $state_settings_ref->{$statename}->{proteinAA} = 1;
        next;
      }

      if(defined $state_settings_ref->{$1}){
        main'warning("State $1 defined more than once in file $cfgfile, tools from first definition may be overriden");
      }
      else{
        $statename = $1;
        $state_settings_ref->{$statename} = {};
      }
    }
    # e.g. spliceorf
    elsif(/^(\S+)\s+(\d+)\s+(\d+)/){
      #next if $1 eq "frame";
      $state_settings_ref->{$statename}->{$1} = [$2, $3];
    }
    elsif(/^(\S+)\s+(\d+)/){
      #next if $1 eq "frame";
      $state_settings_ref->{$statename}->{$1} = $2;
    }
    elsif(/^(graphics)\s+(\S+)/){
      $state_settings_ref->{$statename}->{$1} = $2;
    }
  }
  close(CFG);
  $/ = $tmpln;
}
1;

__END__
