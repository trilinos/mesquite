#!/bin/sh
# A simple script to run configure in the following manner:
# 1. If a config.status exists in the current directory,
#    configure is rerun with the options contained in config.status.
# 2. If config.status doesn't exist and a string with configure
#    options is specified as an argument to reconfig.sh,  
#    configure is executed with these options.
# 3. If config.status doesn't exist and reconfig.sh is run 
#    without arguments, the value of the CCA_CONFIGURE_FLAGS
#    environment variable is used to specify the options to configure.

# One string parameter can be used to specify a string of 
# configure options
if [ -e config.status ]; then
  command=`awk '/^# *\.?\/?configure/  {RS = " "; $1 = ""; \
	   print $0; exit}' config.status ; `;
elif [ "XXX$1" == "XXX" ]; then
  command="./configure $CCA_CONFIGURE_FLAGS";
else 
  command="./configure $1";
fi

echo "Running $command"
$command