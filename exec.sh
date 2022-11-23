#!/bin/bash

name_="ex1"

rm -rf *.o *.gch

g++ -c CRunDec.cpp $name_.cc
g++ -o $name_ $name_.o CRunDec.o

./$name_ > $name_.log

echo "End session."
