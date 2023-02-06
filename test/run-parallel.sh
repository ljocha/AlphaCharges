#!/bin/bash

versions=3,4
phrange=6,8.1,1
par=1


while getopts hv:p:n:s: OPT; do case $OPT in
	v) versions="$OPTARG" ;; 
	p) phrange="$OPTARG" ;;
	n) par=$OPTARG;;
	s) server="-s $OPTARG" ;;
	h|?) echo usage: $0 [-v versions] [-p phrange] [-n parallel] idsfile >&2 ; exit 1;;
esac; done

shift $(($OPTIND - 1))

ids="$1"
lines=$(wc -l "$ids" | cut -d' ' -f1)
chunk=$(($lines / $par))

if [ -f $ids.out ]; then
	echo $ids.out exists, back it up or remove first
	exit 1
fi

line=0
p=0
while read id; do
	if [ $(($line % $chunk)) -eq 0 ]; then
		p=$(($p + 1))
		part=$ids-part$(printf %03d $p)
		rm -f $part
	fi
	echo $id >>$part
	line=$(($line + 1))
done <$ids

dir=$(dirname $0)

trap 'kill $(jobs -p); wait' EXIT
for p in $(seq 1 $par); do
	part=$ids-part$(printf %03d $p)
	$dir/ac-stresstest.py $server -v "$versions" -p "$phrange" $part >$part.out 2>$part.err &
done

wait

trap '' EXIT

rm -f $ids.out
for p in $(seq 1 $par); do
	part=$ids-part$(printf %03d $p)
	cat $part.out >>$ids.out
	rm $part.out
done

