#!/bin/bash

for address in $1;
do

	for i in `find $address -iname result.txt`; do
		echo $i `cat $i | grep '^vertices'` `cat $i | grep 'read time'` `cat $i | grep 'running_time'` `cat $i | grep 'num_rounds'` `cat $i | grep 'total_size'` `cat $i | grep 'total_weight'`;
	done > tmp1.txt;

	cat tmp1.txt | sort | awk -F "[/ ]" '{print $2, $3, $4, $5, $6, $7, $10, $12, $14, $17, $20, $23, $25, $27}' > $address/results.txt;

	sed -i '1s/^/graph algo epsilon k cores repeat v local_e cross_e read_time running_time rounds matching_size matching_weight\n/' $address/results.txt
	
	rm tmp1.txt
 done


