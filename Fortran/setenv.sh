#!/bin/bash
###############################################################################
# Copyright 2020 Ralf Schneider
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################
#
#------------------------------------------------------------------------------
# Setup of the build and execution environment for the 
# HLRS - Material structure process chain
#
# last edited by : Ralf Schneider < schneider@hlrs.de>
#             on : Sep 18th 2019
#
# Prerequisites are :
# --------------------
# + Working mpi installation with 
#   - mpi-compilers
#   - mpirun
# + Working installation of METIS with 64Bit index length
# + Working installation of PETSc with 64Bit index length
#------------------------------------------------------------------------------
#
prefix=$(dirname $BASH_SOURCE)
#
usage ()
{
    echo ""
    echo "Usage:"
    echo "     setenv.sh SYSTEM [silent]"
    echo ""
    echo "Description:"
    echo "     Sets the environment to run the material structure process chain"
    echo "     for SYSTEM."
    echo ""
    echo "Parameters:"
    echo "     SYSTEM : Valid values are 'zeus' or 'vulcan'."
    echo "     silent : Any value for arg 2 supresses all messages"
    echo ""           
}
#
#------------------------------------------------------------------------------
system=$1
#
if [ -z $1 ]; then
    usage
else
    #
    if [ -z $2 ]; then
	echo "============================================================"
	echo "== Setting environment for system : "$1
	echo "=="
    fi
    #
    sys_set=0
    for sys_file in `ls --color=never ${prefix}/auxiliaries/system_environments`
    do
	system=`basename -s .sh $sys_file`
	if [ "$system" == "$1" ]; then
	    source ${prefix}/auxiliaries/system_environments/${sys_file}
	    sys_set=1	    
	fi
    done

    if [ $sys_set == 0 ]; then
       echo ""
       echo "!!! System $1 is currently not supported !!!"
       usage
    else
	#
	# ----------------------------------------
	# struct-process Environment
	#
	# Basepath -------------------------------
	export SP_PREFIX=${prefix}/
	#
	# PATH extensions ------------------------
	export PATH=${prefix}/bin:$PATH
	#
	if [ -z $2 ]; then
	    echo "== Done"
	    echo "============================================================"
	fi
    fi
fi
