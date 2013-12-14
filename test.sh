#!/bin/bash
index=0
while read line; do
    array[$index]="$line"
    index=$(($index+1))
done < names.txt
for ((a=0; a < ${#array[*]}; a++))
do
for ((b=a+1; b < ${#array[*]}; b++))
do
./a.out ${array[$a]} ${array[$b]} 
done
done
	

