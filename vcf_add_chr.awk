#!/bin/awk -f

{if (substr($0,1,1) == "#") print $0; else print "chr"$0}
