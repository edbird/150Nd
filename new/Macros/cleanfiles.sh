#!/bin/bash
FILE="Nd150_2e_P2.root"
if [ -e $FILE ]
then
    rm $FILE
fi
FILE="Nd150_loglikResults.root"
if [ -e $FILE ]
then
    rm $FILE
fi
if [ -e of_numberofeventsafterfit.txt ]
then
    rm of_numberofeventsafterfit.txt
fi
#rm of_numberofeventsafterfit.txt
