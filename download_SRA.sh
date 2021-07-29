#!/usr/bin/env bash

# Created: 06-28-2021
# Updated: 06-28-2021
# 
# Script for downloading a list of SRA files
# EX:
# ./fetch_SRA.sh -o DATA/Weissman_2015_BMDC/ -s DATA/Weissman_2015_BMDC/Weissman_2015_BMDC_SRARunTable.txt

usage() { echo "Usage: $0 [-s <sra_table>] [-o <output_directory>]" 1>&2; exit 0; }

while getopts ":s:o:" arg; do
    case "${arg}" in
        s)
            s=${OPTARG}
            ;;
        o)
            p=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

echo "SRA Table = ${s}"
echo "Output Directory = ${o}"

if [ -z "${s}" ] || [ -z "${o}" ]; then
    usage
fi

sed 1d ${s} | while IFS=, read -r SRA unwanted
do 
	echo "Retrieving $SRA"

	/home/eric/usr_bin/sratoolkit.2.11.0-ubuntu64/bin/fasterq-dump.2.11.0 \
		$SRA \
		--outdir ${o} \
		--temp ${o}/temp \
		--threads 6 \
		--details \
		--progress
done
