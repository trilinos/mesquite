#!/bin/tcsh


echo "This script will add the LGPL license to ALL .hpp and .cpp files in and under this directory. It will NOT check whether the LGPL was already included. Do you want to proceed ? [y/n]"
set proceed = $<

if ($proceed == y) then
foreach i (`find . -regex ".*.hpp"`)
  echo "adding LGPL to $i"
  cat LGPL $i > $i.tmp
  mv -f $i.tmp $i
end

foreach i (`find . -regex ".*.cpp"`)
  echo "adding LGPL to $i"
  cat LGPL $i > $i.tmp
  mv -f $i.tmp $i
end
else
echo "aborting"
endif
