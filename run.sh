#!/bin/bash

./compile.sh

nohup /usr/bin/time -v ./main > log.txt 2>&1 &

disown

exit
