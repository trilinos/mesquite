#!/bin/tcsh


echo "Do you want to remove all CVS directories in and under this directory ? [y/n]" 
set proceed = $<


if ($proceed == y) then
foreach i (`find . -regex ".*/CVS"`)
  echo "removing directory $i"
  rm -rf $i
end
else
echo "aborting"
endif
