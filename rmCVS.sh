#!/bin/tcsh

foreach i (`find . -regex ".*/CVS"`)
  echo "removing directory $i"
  rm -rf $i
end
