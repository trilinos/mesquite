######################################################################################
# Extract the value of a variable from a makefile fragment.
# Arguments:
#   - The path of the makefile fragment
#   - The name of the makefile variable
#   - Action on success (the shell variable "make_val" will contain the value
#         of the variable)
#   - Action on failure
#######################################################################################
AC_DEFUN([SNL_MAKE_INC_VAR], [
make_val=
snl_makefile="snl_check.mak"
snl_resultfile="snl_check.out"
rm -f $snl_makefile

if test ! -f $1 ; then
  AC_MSG_WARN([File not found: $1])
  $4
else
cat >$snl_makefile <<SNL_END_OF_MAKEFILE
default:
	@echo "\$($2)" > $snl_resultfile

include $1
SNL_END_OF_MAKEFILE
if make -f $snl_makefile > /dev/null 2>&1; then
  make_val=`cat $snl_resultfile`
  rm -f $snl_makefile
  rm -f $snl_resultfile
  $3
else
  rm -f $snl_makefile
  rm -f $snl_resultfile
  $4
fi
fi
])

##############################################################################
# Optinally configure support for an ITAPS API
#
# Arguments:
#   1) The interface name in lowercase (e.g. imesh)
#   2) The interface name in uppercase (e.g. IMESH)
#   3) The interface name in preferred case (e.g. iMesh)
#   4) Default value (yes/no)
#
# Sets and declares substitutions for:
#    ENABLE_IFACE=yes/no
#    WITH_IFACE_IMPL=yes/no
#    IFACE_DEFS=empty/path to defs file
#    IFACE_LIBS=empty/ld args for linking implementation
# Sets AUTOMAKE conditionals for ENABLE_IFACE and WITH_IFACE_IMPL
##############################################################################

AC_DEFUN([ITAPS_API], [

AC_ARG_ENABLE( [$1], 
[AC_HELP_STRING([--disable-$1],
  [Do not build support for ITAPS $1 interface.])],
AC_HELP_STRING([--enable-$1=DIR],
  [Build support for ITAPS $1 interface.  Optionally specify implementation.])
  [ENABLE_$2=$enableval],[ENABLE_$2=yes])

WITH_$2_IMPL=no
$2_DEFS=
$2_LIBS=
if test "xno" != "x$ENABLE_$2"; then
  if test "xyes" != "x$ENABLE_$2"; then
    for subdir in . lib include; do
      if test -f "$ENABLE_$2/$subdir/$3-Defs.inc"; then
        $2_DEFS=$ENABLE_$2/$subdir/$3-Defs.inc
        break
      fi
    done
    if test "x" = "x${$2_DEFS}"; then
      AC_MSG_ERROR("$3-Defs.inc not found in $ENABLE_$2")
    else
      SNL_MAKE_INC_VAR( [${$2_DEFS}], [$2_LIBS], [$2_LIBS="$make_val"] )
      WITH_$2_IMPL=yes
    fi
    ENABLE_$2=yes
  fi
fi

AC_SUBST(ENABLE_$2)
AC_SUBST(WITH_$2_IMPL)
AC_SUBST($2_DEFS)
AC_SUBST($2_LIBS)
AM_CONDITIONAL([ENABLE_$2],[test "xyes" = "x$ENABLE_$2"])
AM_CONDITIONAL([WITH_$2_IMPL],[test "xyes" = "xWITH_$2_IMPL"])

])

