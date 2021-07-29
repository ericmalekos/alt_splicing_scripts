#!/usr/bin/env python3

'''
Utility for matching genomic coordinates - from BED file - to the closest annotated genes.

Requires a BED file and sorted GTF file. To sort GTF use:

    bedtools sort -i <gtf.in> | awk '$3 == "gene" {print}' > <genes_only_sorted.gtf>

TODO:


'''

import pybedtools as pbt
import argparse
import sys

def parseCommandLine():
    parser = argparse.ArgumentParser(description='Convert alternative splicing results to BED format')
    parser.add_argument('-p', '--program', type=str, help='Splicing Program', choices=['rmats', 'drimseq', 'juncbase'])
    parser.add_argument('-b', '--bed', type=str, help='File Name')
    parser.add_argument('-g', '--gtf', type=str, help='File Name')
    parser.add_argument('-s', '--stranded', help='Ignore novel event types in DRIMseq/juncBASE', default=False, action='store_true')
    return parser.parse_args()

def fillBlanksDRIMSeq(bedIn, sortedGTFIn, stranded=False, startPattern="chr"):
    ''' Some juncBASE/DRIMSeq entries had no geneID i.e. no "XLOC...". After converting to BED
    with SpliceResultsToBed, the DRIMSeq files that have no geneID instead have "chr:##-##;[K,N]"
    in the 4th column position.
    This method works by identifying those missing geneIDs and replacing them with whatever geneIDs
    either overlap or are closest to them in the GTF file.
    '''

    # TODO This gave some strange results: extra numbers and some missing genes

    genes = pbt.BedTool(sortedGTFIn)
    with open(bedIn, 'r') as bed:
        for line in bed.readlines():
            split_line = line.strip().split("\t")
            if split_line[3].startswith(startPattern):
                bed_entry = pbt.BedTool(line, from_string=True)
                closest = bed_entry.closest(genes, s=stranded, d=True)
                overlapping_genes = []
                for entry in closest:
                    overlapping_genes.append(str(entry).split("gene_name")[1].split(";")[0].strip().strip("\""))
                split_line[3] = ",".join(overlapping_genes)
                print("\t".join(split_line))
            else:
                print(line.strip())
def main():

    args = parseCommandLine()
    # if args.program == "rmats":
    #     rMATStoBED(preprocessRMATS(fileIn=args.input, FDR=args.padj, inclevel=args.inclusionLevel, 
    #                     readSupport=args.minCount, save=args.saveAs), eventType=args.event)

    if args.program == "drimseq":
        fillBlanksDRIMSeq(bedIn=args.bed, sortedGTFIn=args.gtf, stranded=args.stranded)
        
    # if args.program == "juncbase":
    #     JuncBASEtoBED(preprocessJuncBASE(fileIn=args.input, adjP=args.padj, save=args.saveAs, k=args.known, eventType=args.event))


if __name__ == "__main__":
    main()