#!/bin/bash
#command="./joblaunch.sh"
evalcmd="eval $(~/snemo-sdk/bin/brew shellenv)"
command="root"
arguments="-b ../newLogLikFitter.C"
numcores=8
for i in {0..301}
do

    sleep 1
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
