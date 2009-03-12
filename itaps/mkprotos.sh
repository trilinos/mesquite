#!/bin/sh

PFX="$1"
INPUT="$2"
OUTPUT="$3"

if test "x" = "x$SED"; then
  SED=`which sed`
fi

EXPR="s/^[[:space:]]*void[[:space:]][[:space:]]*${PFX}_\([a-z][_a-zA-Z0-9]*\)[[:space:]]*(.*\$/${PFX}_\1/p"

echo '#include "iBase_FCDefs.h"' > $OUTPUT
echo '#ifdef FC_FUNC_' >> $OUTPUT
echo >> $OUTPUT
for func in `$SED -n "$EXPR" $INPUT`; do
  lower=`echo $func | tr '[:upper:]' '[:lower:]'`
  upper=`echo $func | tr '[:lower:]' '[:upper:]'`
  echo "#define $func FC_FUNC_( $lower, $upper )" >> $OUTPUT
done
echo >> $OUTPUT
echo "#endif" >> $OUTPUT
