#!/bin/bash

rm -f signalserver
rm -f signalserverHD
g++ -Wall -O3 -s -lm -fomit-frame-pointer itwom3.0.cpp models.cpp main.cpp -o signalserver
./signalserver
