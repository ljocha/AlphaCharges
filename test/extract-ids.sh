#!/bin/bash 

# requires google-cloud-cli, see https://cloud.google.com/storage/docs/gsutil_install

# run "gcloud auth login" first

for i in $(seq 1 22); do 
	n=$(printf %03d $i)
	sed 's/^AF-\([^-]*\)-.*$/\1/' manifest-model_v4_cif-part-$n.csv >ids-part-$n.txt
done
