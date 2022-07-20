use Conf;

sub main'fatal{
  die "$0 FATAL: ", join("", @_), "\n";
}

sub main'warning{
  warn "$0 WARNING: ", join("", @_), "\n";
}

sub main'notify{
  my ($project, @message, $project_settings_ref) = @_;

  if(not defined $project_settings_ref){
    my %ps = ();
    $project_settings_ref = \%ps;
    &Config::project($project, $project_settings_ref)
  }

  my $mailprog = "/bin/mail";
  if(!open(MAIL, "| $mailprog ". $$project_settings_ref{'notify'})){
    warning("Couldn't mail message to ", $$project_settings_ref{'notify'},
	    "using $mailprog: $!\nMessage content:\n", @message);
    return;
  }

  print MAIL "From: Magpie project $project\n",
             "Subject: Automated notification\n\n", @message;
  
  close(MAIL);
}

1;

__END__

