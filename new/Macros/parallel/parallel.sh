#!/bin/bash
#command="./joblaunch.sh"
evalcmd="eval $(~/snemo-sdk/bin/brew shellenv)"
command="root"
arguments="-b ../newLogLikFitter.C"
numcores=11
for i in {0..12}
do

    sleep 10
    while [ $( pgrep -c -x "$command" ) -ge "$numcores" ]
    do
        kill -SIGSTOP $$
    done

    (
        eval $(~/snemo-sdk/bin/brew shellenv)
        root -b "newLogLikFitter.C($i)"
        kill -SIGCONT $$
    ) &

done
