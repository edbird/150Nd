#!/bin/bash
#command="./joblaunch.sh"
#evalcmd="eval $(~/snemo-sdk/bin/brew shellenv)"
command="root"
arguments="-b ../newLogLikFitter.C"
numcores=12
i=1 # parallel jobs start from 1
eval $(~/snemo-sdk/bin/brew shellenv)
while [[ $i -lt $numcores ]]
do
    echo "EXEC $i"
    root -b "newLogLikFitter.C($i)"
    ((i = i + 1))
done
#for i in {0..12}
#do
#
#    sleep 10
#    while [ $( pgrep -c -x "$command" ) -ge "$numcores" ]
#    do
#        kill -SIGSTOP $$
#    done
#
#    (
#        eval $(~/snemo-sdk/bin/brew shellenv)
#        root -b "newLogLikFitter.C($i)"
#        kill -SIGCONT $$
#    ) &
#
#done
