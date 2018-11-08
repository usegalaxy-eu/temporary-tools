#!/usr/bin/env bash
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#

export PATH=$PATH:$(dirname $0)

toolExample.pl $*

exit 0
