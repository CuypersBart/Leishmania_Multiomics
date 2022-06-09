import os
import numpy as np
import pandas as pd
import math
from scipy.stats import gmean, mannwhitneyu, ttest_ind
from IPython.display import IFrame
import random

def checkSNP(strain1, strain2, gene, df):
    output = []
    variants = df.loc[df['Genename'] == gene]
    for index, row in variants.iterrows():
        if row[strain1] != row[strain2]:
            output.append([row[strain1], row[strain2], row['SNP_type'], row['SNP_impact'], row['SNP_substitution_base'], row['SNP_substitution_amino'], row['SNP_pos_transcript'], row['SNP_pos_protein']])
    return output

def averager(q, mindetect, mean = 'arithmetic'):
    y = [float(i) for i in q]
    z = []
    for i in y:
        if float(i) != float(0):
            z.append(i)
    if len(z) >= mindetect:
        if mean == 'GM':
            return float(gmean(z))
        else:
            return float(np.average(z))
    else:
        return float("NaN")

def readomic(omicfile, *args):
    lines = -1
    for ar in args:
        lines = ar
    with open (omicfile) as myfile:
        header = myfile.readline()
        samples = header.strip('\n\r').replace('"', '').split('\t')[1:]
        mydict = {}
        while lines != 0:
            i = myfile.readline()
            if not i:
                break
            i = i.strip('\n\r').replace('"', '').split('\t')
            mydict[str(i[0])] = {}
            for j in range(1, len(samples) + 1):
                mydict[str(i[0])][samples[j - 1]] = i[j]
            lines = lines - 1
    return [mydict, samples]

def dicttopandas(df):
    return pd.DataFrame.from_dict(df, orient='index', dtype = float)

def normalise(df):
    #add pseudocounts for remaining zeros
    df = df.replace(0, 1)
    #take 2 log
    df = np.log2(df)
    #divide column by column mean
    for column in df:
        df[column] = df[column] - np.average(df[column])
    return(df)

def genefoldchange(strain1, strain2, pd, datatype, chr = 'all', mindect = 3):
    normexp = {}
    strains = [strain1, strain2]
    for gene, row in pd.iterrows():
        values = [[],[]]
        proceed = False
        currentchrom = gene.split('_')[1][0:2]
        if chr == 'all':
            proceed = True
        elif currentchrom == chr[2:]:
            proceed = True
        if proceed == True:
            for i in [0,1]:
                if '275' not in strains[i]:
                    Sample = strains[i]
                    values[i] = averager([float(row[Sample + '.1']), float(row[Sample + '.2']),
                                          float(row[Sample + '.3']), float(row[Sample + '.4'])], mindect)
                else:
                    Sample = strains[i]
                    values[i] = averager([float(row[Sample + '.1']), float(row[Sample + '.2']),
                                          float(row[Sample + '.3'])], mindect)
            if datatype == "LOG":
                normexp[gene] = values[1] - values[0]
            else:
                normexp[gene] = np.log2(values[1] / values[0])
    return normexp

def foldchangePT(Groups, diploidref, pd, datatype, select, randomref = False):
    normexp = {}
    for gene, row in pd.iterrows():
        chr = 'Ld' + gene.split('_')[1][0:2]
        if chr in diploidref:
            if chr not in normexp:
                normexp[chr] = {}
            normexp[chr][gene] = {}
            refvalues = []
            tempsamples = []
            for refsample in diploidref[chr]:
                if float(row[refsample[0]]) > 0.0:
                    for refsample in diploidref[chr]:
                        refvalues.append(row[refsample[0]] / refsample[1])
                        tempsamples.append(refsample)
            if refvalues == []:
                refvalue = 0
                out = "NotinHere"
            else:
                if randomref == True:
                    refvalue = random.choice(refvalues)
                    out = tempsamples[refvalues.index(refvalue)][0].split('.')[0]
                else:
                    out = 'NA'
                    refvalue = np.median(refvalues)
            for Sample in Groups:
                if out not in Sample:
                    if select == 'all':
                        averageexp = averager([float(row[Sample + '.1']), float(row[Sample + '.2']),
                                              float(row[Sample + '.3']), float(row[Sample + '.4'])], 2, 'GM')
                    else:
                        averageexp = float(row[Sample + select])
                    if datatype == "LOG":
                        normexp[chr][gene][Sample] = averageexp - refvalue
                    else:
                        if refvalue != 0 and averageexp != 0:
                            normexp[chr][gene][Sample] = np.log2(averageexp / refvalue)
                        else:
                            normexp[chr][gene][Sample] = float("NaN")
                else:
                    normexp[chr][gene][Sample] = float("NaN")
    return normexp

def foldchange(Groups, diploidref, pd, datatype, select, randomref = False):
    normexp = {}
    for gene, row in pd.iterrows():
        chr = 'Ld' + gene.split('_')[1][0:2]
        if chr in diploidref and row[0]:
            if chr not in normexp:
                normexp[chr] = {}
            normexp[chr][gene] = {}
            refvalues = []
            tempsamples = []
            for refsample in diploidref[chr]:
                tempref = row[refsample[0]]
                tempref = 2 ** tempref
                tempref = tempref / refsample[1]
                tempref = np.log2(tempref)
                refvalues.append(tempref)
                tempsamples.append(refsample)
            if tempsamples == []:
                for Sample in Groups:
                    normexp[chr][gene][Sample] = float("NaN")
            else:
                if randomref == True:
                    refvalue = random.choice(refvalues)
                    out = tempsamples[refvalues.index(refvalue)][0].split('.')[0]
                else:
                    refvalue = np.median(refvalues)
                    out = 'NA'
                for Sample in Groups:
                    if out not in Sample:
                        if select == 'all':
                            if '275' not in Sample:
                                averageexp = averager([float(row[Sample + '.1']), float(row[Sample + '.2']),
                                                      float(row[Sample + '.3']), float(row[Sample + '.4'])], 4)
                            else:
                                averageexp = averager([float(row[Sample + '.1']), float(row[Sample + '.2']),
                                                      float(row[Sample + '.3'])], 3)
                        else:
                            averageexp = float(row[Sample + select])
                        if datatype == "LOG":
                            normexp[chr][gene][Sample] = averageexp - refvalue
                        else:
                            normexp[chr][gene][Sample] = np.log2(averageexp / refvalue)
                    else:
                        normexp[chr][gene][Sample] = float("NaN")
    return normexp

def chromaverager(normexp, groups):
    chromav = {}
    for chr in normexp:
        chromav[chr]={}
        for sample in groups:
            chromav[chr][sample] = []
        for gene in normexp[chr]:
            for sample in normexp[chr][gene]:
                if math.isnan(normexp[chr][gene][sample]) == False and math.isinf(normexp[chr][gene][sample]) == False:
                    chromav[chr][sample].append(normexp[chr][gene][sample])
        for sample in chromav[chr]:
            chromav[chr][sample] = [np.median(chromav[chr][sample]), len(chromav[chr][sample])]
    return(chromav)

def ploidyperchrom(ploidydict):
    PloidyFCperchrom = {}
    for chr in ploidydict:
        PloidyFCperchrom[chr] = {}
        for sample in ploidydict[chr]:
            if '.1' in sample:
                PloidyFCperchrom[chr][sample.replace('.1','')] = np.log2(float(ploidydict[chr][sample]))
    return PloidyFCperchrom

def transformer(mylist):
    return("\t".join([str(2**(x)) for x in mylist]))