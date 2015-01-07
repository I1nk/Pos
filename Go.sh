#!/bin/bash

clear
gcc -c -Wall Pos.c -O2
gcc -static Pos.o -lm -o Pos.out -O2


