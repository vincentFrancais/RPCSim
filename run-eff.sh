#!/bin/bash

let "i=0"
let "hv=51"

while [ $i -lt 2 ]
do
	time ./sim config/calice.xml $hv
	let "hv += 1"
	let "i += 1"
done 
	
