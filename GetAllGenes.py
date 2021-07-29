#!/usr/bin/env python3

'''

TODO:


'''

import argparse
from collections import defaultdict
import xlsxwriter

rMATSEventDict = {"Alternative 3' Splice Site selection from Junction Count": "A3SS_JC",
                "Alternative 3' Splice Site selection from Junction and Exon Count": "A3SS_JCEC",
                "Alternative 5' Splice Site selection from Junction Count": "A5SS_JC",
                "Alternative 5' Splice Site selection from Junction Count and Exon Count": "A5SS_JCEC",
                "Mutually Exclusive Exon from Junction Count": "MXE_JC",
                "Mutually Exclusive Exon from Junction Count and Exon Count": "MXE_JCEC",
                "Retained Intron from Junction Count": "RI_JC",
                "Retained Intron from Junction Count and Exon Count": "RI_JCEC",
                "Skipped Exon from Junction Count" : "SE_JC",
                "Skipped Exon from Junction Count and Exon Count" : "SE_JCEC"}

DRIMSeqEventDict = {'alternative_donor': 'A5SS', 'alternative_acceptor': 'A3SS',
                    'cassette': 'SE', 'intron_retention': 'RI', 'mutually_exclusive': 'MXE',
                    'coord_cassette': 'COORD', 'alternative_first_exon': 'AFE', 'alternative_last_exon': 'ALE'}

def parseCommandLine():
    parser = argparse.ArgumentParser(description='Convert alternative splicing results to BED format')
    parser.add_argument('-r', '--rmats', type=str, help='Space separated list of rMATS files', nargs='*')
    parser.add_argument('-d', '--drimbed', type=str, help='Space separated list of BED files from DRIMSeq', nargs='*')
    parser.add_argument('-c', '--combinedrMATS', type=str, help='Use this option for combined rMATS files')
    parser.add_argument('-o', '--output', type=str, help='Prefix of output XLSX file')

    return parser.parse_args()

def namesFromCombinedRMATS(rmatsIn, sep=","):
    rmats_dict = {}
    cur_key = ""
    with open(rmatsIn, 'r') as rmats:
        for entry in rmats.readlines():
            if entry.startswith("#") and len(entry.strip()) == 1:
                continue
            elif entry.startswith("#") and len(entry.strip()) > 1:
                cur_key = rMATSEventDict[entry.split("\t")[-1].strip()]
                rmats_dict[cur_key] = []
            else:
                genes = entry.split("\t")[2].split(sep)
                for gene in genes:
                    if gene not in rmats_dict[cur_key] and gene != "GeneID":
                        rmats_dict[cur_key].append(gene)
                    rmats_dict[cur_key].sort()
    
    return rmats_dict

def nameFromDrimBED(bedIn, sep=","):
    key = bedIn.split("/")[-1].split(".")[0][12:]
    drimBedDict={}
    drimBedDict[key] = []
    with open(bedIn, 'r') as bed:
        for entry in bed.readlines():
            names = entry.split("\t")[3].split(sep)
            for name in names:
                if name not in drimBedDict[key]:
                    drimBedDict[key].append(name)
                drimBedDict[key].sort()

    return drimBedDict

def combinedValues(dictIn):
    all_values = []
    for value in dictIn.values():
        for v in value:
            if v not in all_values: all_values.append(v)
    return sorted(all_values)

def writeDictToXLSX(output, dictin):
    workbook  = xlsxwriter.Workbook(output + ".xlsx")
    for tab, values in dictin.items():
        worksheet = workbook.add_worksheet(tab)

        for row_ndx, value in enumerate(values):
            worksheet.write(row_ndx, 0, value)

    workbook.close()

def main():

    args = parseCommandLine()
    combinedDict={}
    if args.drimbed:
        for drim in args.drimbed:
            combinedDict.update(nameFromDrimBED(drim))
        combinedDict["All DRIMSeq"] = combinedValues(combinedDict)

    #print(combinedDict)
    if args.combinedrMATS:
        rMATSDict = namesFromCombinedRMATS(args.combinedrMATS)
        rMATSDict["All rMATS"] = combinedValues(rMATSDict)
        combinedDict.update(rMATSDict)

    combinedDict["rMATS-DRIMSeq Union"] = combinedValues(combinedDict)

    writeDictToXLSX(output=args.output, dictin=combinedDict)


if __name__ == "__main__":
    main()
