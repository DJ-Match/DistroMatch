for i in $1 ; \
do \
	echo `find $i -iname gen_map.sh | wc -l` \
		`find $i -iname result.txt | wc -l` $i ;  \
done
