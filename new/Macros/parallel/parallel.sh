#!/bin/bash
#command="./joblaunch.sh"
#evalcmd="eval $(~/snemo-sdk/bin/brew shellenv)"
command="root"
arguments="-b ../newLogLikFitter.C"
numcores=12
i=1 # parallel jobs start from 1
#eval $(~/snemo-sdk/bin/brew shellenv)
while [[ $i -le $numcores ]]
do
    echo "EXEC $i"
    (
        eval $(~/snemo-sdk/bin/brew shellenv)
        root -b "newLogLikFitter.C($i)" > "cout_$i.txt"
        #kill -SIGCONT $$
    ) &
    ((i = i + 1))
    sleep 5
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
