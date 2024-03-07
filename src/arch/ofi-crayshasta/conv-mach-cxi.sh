
# For libfabric
#If the user doesn't pass --basedir, use defaults for libfabric headers and library
if test -z "$USER_OPTS_LD"
then
    if test -z "$LIBFABRIC"
    then
	CMK_INCDIR="$CMK_INCDIR -I/usr/include/"
	CMK_LIBDIR="$CMK_LIBDIR -L/usr/lib64/"
    else
	CMK_INCDIR="$CMK_INCDIR -I$LIBFABRIC/include/"
	CMK_LIBDIR="$CMK_LIBDIR -L$LIBFABRIC/lib64/"
    fi
fi

# For cray-pmi
if test -n "$CRAY_PMI_PREFIX"
then
    CMK_INCDIR="$CMK_INCDIR -I$CRAY_PMI_PREFIX/include"
    CMK_LIBDIR="$CMK_LIBDIR -L$CRAY_PMI_PREFIX/lib"
fi

CMK_LIBS="$CMK_LIBS -lfabric"
# Use PMI2 by default on Cray systems with cray-pmi
. $CHARMINC/conv-mach-slurmpmi2.sh

# For runtime
CMK_INCDIR="$CMK_INCDIR -I./proc_management/"
