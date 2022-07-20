#!/usr/local/bin/perl

use GDBM_File;
tie(%db, "GDBM_File", $ARGV[0], O_RDONLY, 0644) or die "Cannot open database: $!\n";
#print join("\n", keys %db), "\n";
  print $db{$ARGV[1]}, "\n";
