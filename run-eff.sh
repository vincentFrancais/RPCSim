#!/bin/bash

let "i=0"
let "hv=40"

while [ $i -lt 20 ]
do
	time ./sim $hv $hv
	let "hv += 1"
	let "i += 1"
done 
	
