
package Tool;

use strict;
use Constants;
use vars qw($AUTOLOAD);

sub new {
  shift;
  my $self = bless {};
  my $project = shift;
  $self->{'name'} = shift;
  if($_[0] != $Constants::protein && $_[0] != $Constants::dna &&
     $_[0] != $Constants::hmm && $_[0] != $Constants::msa){
    warn "Third argument must either Constants::protein, Constants::dna, Constants::msa or Constants::hmm was $_[0]\n";
    return undef;
  }
  $self->{'input_molecule'} = shift;

  if($_[0] !~ /^(local|\S+\@\S+|http)$/i){
    warn "Fourth argument must be one of \"local\", an e-mail or HTTP address, was $_[0]\n";
    return undef;
  }
  $self->{'method'} = shift;

  if($_[0] != $Constants::facttool && $_[0] != $Constants::texttool){
    warn "Fifth argument must be either Constants::facttool of Constants::texttool\n"; 
    return undef;
  }
  $self->{'evid_type'} = shift;

  if($_[0] != $Constants::dna and $_[0] != $Constants::protein and
     $_[0] != $Constants::motif and $_[0] != $Constants::orf and $_[0] != $Constants::model and
     $_[0] != $Constants::trna and $_[0] != $Constants::repeat){
    warn "Sixth argument must be one of Constants::(dna|protein|motif|model|orf), was $_[0]\n";
    return undef;
  }
  $self->{'sim_type'} = shift;

  if($_[0] != $Constants::allscope && $_[0] != $Constants::interorfscope &&
     $_[0] != $Constants::orfscope){
    warn "Seventh argument must be one of Constants::(all|interorf|orf)scope, was $_[0]\n";
    return undef;
  }
  $self->{'scope'} = shift; 

  $self->{'priority'} = shift;

  $self->{'run_priority'} = shift;

  $self->{'run_hosts'} = shift;

  # How many sequences can we send into the tool at once?
  $self->{'run_chunks'} = shift;

  if($self->{'evid_type'} == $Constants::facttool){
    if(not $_[0]){
        warn "Missing nineth argument, parser name, required when tool produces facts";
        return undef;
    }
    $self->{'parser'} = shift;
    if($self->{'parser'} =~ /,/){
      my @parsers = split /,/, $self->{'parser'};
      my %parsers;
      my $spec;
      while($spec = pop @parsers){
        if($spec =~ /^(\S+):(\S+)$/){
          $parsers{$1} = $2;
        }
        else{
          warn "Could not parse tool parser spec (skipping): $spec\n";
        }
      }
     
      $self->{'parser'} = \%parsers;

     
    }
  }

  $self->{'format'} = shift;

  $self->{'_permitted'} = { project => 1, name => 1, input_molecule => 1,
                            method => 1, evid_type => 1, sim_type => 1, run_chunks => 1,
                            run_priority => 1, run_hosts => 1, scope => 1, priority => 1, 
                            parser => 1, format => 1};
  return $self;
}

sub AUTOLOAD{
  my $self = shift;
  my $name = $AUTOLOAD;
  $name =~ s/^.*://;
  return if $name eq "DESTROY";
  die "Can't access '$name' field in object ", ref $self
    unless exists $self->{_permitted}->{$name};

  if(@_){
    return $self->{$name} = shift;
  }
  else{
    return $self->{$name};
  }
}

1;
__END__

