#!/usr/local/bin/perl

$config_flags = $ENV{"CCA_CONFIGURE_FLAGS"};

$command = "./configure $config_flags";
if (-e "config.status") {
  $command = "echo \"Couldn't parse config.status \"";
  @head = `head config.status`;
  foreach (@head) {
    if (/^#\s(configure\s.*$)/ || /^#\s(\.\/configure\s.*$)/)  {
       $command = "$1";
    }
  }
}
print "Running $command ...\n";
system($command);

