#!/bin/bash 

# requires google-cloud-cli, see https://cloud.google.com/storage/docs/gsutil_install

# run "gcloud auth login" first

for i in $(seq 1 22); do 
	n=$(printf %03d $i)
	gsutil -m cp "gs://public-datasets-deepmind-alphafold-v4/manifests/manifest-model_v4_cif-part-$n.csv" .
done
