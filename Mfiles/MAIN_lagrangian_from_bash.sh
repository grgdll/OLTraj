#!/usr/bin/octave-cli

# run in parallel using: echo {1998..2004..1} | xargs -P7 -n1 echo MAIN_lagrangian_from_bash.sh

# Arg 1 is the year to process
args = argv();
for imonth=2:2
   MAIN_lagrangian_diags(str2num(args{1}),imonth);
end
