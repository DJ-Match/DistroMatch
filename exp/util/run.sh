for address in $1;
do
	touch -a ./${address}/files.txt
	find ./$address -iname gen_map.sh > ./${address}/files.txt;
  	cat ./${address}/files.txt;
done | parallel -j$2 sh;
