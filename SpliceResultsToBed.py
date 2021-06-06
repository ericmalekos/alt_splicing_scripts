
'''
Utility for filtering the statistical analysis output from rMATS and DRIMSeq.
The BED data is sent to stdout, other messages to stderr.

To obtain correctly sorted BED file run as:

    python SpliceResultsToBed.py -p rmats -e RI -i RI.MATS.JCEC.txt | sort -k1,1V -k2,2n -k3,3n

    python SpliceResultsToBed.py -p drimseq -e RI -i drimResults_intron_retention.tsv --known | sort -k1,1V -k2,2n -k3,3n

TODO:
    - implement juncbase option
    - add other event types: currently only intron retention is included
'''


import pandas as pd
import argparse
import sys

from pandas.core import indexing

def parseCommandLine():
    parser = argparse.ArgumentParser(description='Convert alternative splicing results to BED format')
    parser.add_argument('-p', '--program', type=str, help='Splicing Program', choices=['rmats', 'drimseq', 'juncbase'])
    parser.add_argument('-e', '--event', type=str, help='Splicing Event', choices=['A5SS', 'A3SS', 'SE', 'RI', 'MXE', 'COORD', 'AFE', 'ALE'])
    parser.add_argument('-i', '--input', type=str, help='File Name')
    parser.add_argument('-q', '--padj', type=float, help='padj/FDR Cutoff',default=0.1)
    parser.add_argument('-l', '--inclusionLevel', type=float, help='PSI Cutoff for rMATS', default=0.1)
    parser.add_argument('-m', '--minCount', type=int, help='Minimum Read Support for rMATS', default=10)
    parser.add_argument('-s', '--saveAs', type=str, help='Save processed filtered TSV results', default="")
    parser.add_argument('-k', '--known', help='Ignore novel event types in DRIMseq/juncBASE', default=False, action='store_true')
    return parser.parse_args()

def sumReadCounts(entry, sep=","):
    '''Helper function for preprocessRMATS, splits replicate read counts and returns sum'''
    return sum([int(e) for e in entry.split(sep)])

def knownOnly(entry):
    '''Helper function for preprocessDRIMSeq, returns true if the event type reflects a known feature, false if novel'''
    if ';K;' in entry:
        return True
    elif ';N;' in entry:
        return False
    else:
        raise ValueError('Could not determine if feature is Known or Unknown. Investigate')

def preprocessRMATS(fileIn, FDR, inclevel, readSupport, save=""):
    sys.stderr.write("\nProcessing rMATS file: " + fileIn + "\n")
    df = pd.read_csv(fileIn, sep='\t', header=0).sort_values(by='FDR')
    sys.stderr.write("Number of event entries before processing: " + str(df.shape[0]) + "\n")
    df = df[(df.FDR <= FDR) & (abs(df.IncLevelDifference) >= inclevel)]
    df[df.IJC_SAMPLE_1.apply(sumReadCounts) + df.IJC_SAMPLE_2.apply(sumReadCounts) + 
    df.SJC_SAMPLE_1.apply(sumReadCounts) + df.SJC_SAMPLE_2.apply(sumReadCounts)>= readSupport]
    sys.stderr.write("Number of event entries after processing: " + str(df.shape[0])+"\n\n")
    
    #print(df)
    if len(save) > 0:
        df.to_csv(save, sep='\t', index=False)

    return df

def tabSplitRow(row):
    '''Helper function for rMATStoBED, prints tab separated string of row values'''
    print('\t'.join(map(str,row.values)))

def rMATStoBED(df, eventType):

    if eventType == 'RI':
        df[["chr", "upstreamEE", "downstreamES", "GeneID", "FDR", "strand"]].apply(tabSplitRow,axis=1)
    


def preprocessDRIMSeq(fileIn, adjP, save="", k=False):
    sys.stderr.write("Processing DRIMSeq file: " + fileIn + "\n")
    df = pd.read_csv(fileIn, sep='\t', names=[ ' ', 'Event', 'lr', 'df', 'pvalue', 'adj_pvalue' ], skiprows=1).sort_values(by='adj_pvalue')
    sys.stderr.write("Number of event entries before processing: " + str(df.shape[0]) + "\n")
    df = df[(df.adj_pvalue <= adjP)]
    if k:
        sys.stderr.write("Removing Novel Events\n")
        df = df[df.Event.apply(knownOnly)]
    sys.stderr.write("Number of event entries before processing: " + str(df.shape[0]) + "\n")
    
    #print(df)
    if len(save) > 0:
        df.to_csv(save, sep='\t', index=False)
    
    return df

def DRIMtoBED(df):
    events = df['Event'].tolist()

    bed_list = []
    for event in events:
        l = event.split(";")[1].split(":")
        bed_entry = "\t".join([l[0],l[1].split("-")[0],l[1].split("-")[1], event.split(";")[-2]+ ";" + event.split(";")[-1]])
        print(bed_entry)
    

def main():

    args = parseCommandLine()
    if args.program == "rmats":
        rMATStoBED(preprocessRMATS(fileIn=args.input, FDR=args.padj, inclevel=args.inclusionLevel, 
                        readSupport=args.minCount, save=args.saveAs), eventType=args.event)


    if args.program == "drimseq":
        DRIMtoBED(preprocessDRIMSeq(fileIn=args.input, adjP=args.padj, save=args.saveAs, k=args.known))


if __name__ == "__main__":
    main()