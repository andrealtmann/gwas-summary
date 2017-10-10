#! /usr/bin/python

import os, sys, string, gzip

#script produces two files to be used together with PLINK:
#1) one file containg SNP and p-value
#2) second file containing SNP allele and \beta (aka score)

#This script has to:
# 1) identify SNP (rsID) from source file,
#    if the RS number is not provided: convert chr pos info into RS
#    by looking up in a correspondong datbase
# 2) identify the effect allele
# 3) identify the score


def processSumStat(fname, ofbase, varmap):

    tok = fname.strip().split('.')
    if tok[-1] == "gz":
        f = gzip.open(fname)
    else:
        f = open(fname)

    header = f.readline()
    htok   = header.strip().split()
    #get column numbers
    varid = {}
    for v in varmap:
        cn = varmap[v]
        try:
            myid = htok.index(cn)
        except ValueError:
            myid = -1
            sys.stderr.write("Column header '" + v + "' not found. Aborting.\n")
            sys.exit(myid)
        varid[v] = myid
    #print varid
    #we made it till here, let's open the output files
    fpv = open(ofbase + ".pvalue", 'w')
    fpv.write("SNP\tP\n")
    fsc = open(ofbase + ".scores", 'w')
    line = f.readline()
    while line != "":
        tok = line.strip().split()
        retall=False
        try:
            snpname = tok[varid["SNP"]]
            pvalue  = float(tok[varid["pvalue"]])
            vscore  = float(tok[varid["score"]])
            effal   = tok[varid["EA"]]
            retall  = True
        except IndexError:
            sys.stdout.write("Failed to retreive information for SNP: " + snpname + "\n")

        if retall:
            fpv.write(snpname + "\t" + str(pvalue) + "\n")
            fsc.write(snpname + "\t" + effal + "\t" + str(vscore) + "\n")
        line = f.readline()

    f.close()
    fpv.close()
    fsc.close()

def help(sname, varmap):
    sys.stderr.write("USAGE: " + sname + " [options] <infile> <outbase>\n\n")
    sys.stderr.write("options:\n")
    sys.stderr.write("--ea <name>\tcolumn name containg the effect allele, default: " + varmap["EA"] + "\n")
    sys.stderr.write("--pv <name>\tcolumn name containg the p-value, default: " + varmap["pvalue"] + "\n")
    sys.stderr.write("--score <name>\tcolumn name containg the score, default: " + varmap["score"] + "\n")
    sys.stderr.write("--snpname <name>\tcolumn name containg the SNP, default: " + varmap["SNP"] + "\n")
    sys.exit(-1)

### main program ###

varmap = {}
varmap["EA"] = "Effect_allele"
varmap["pvalue"] = "P"
varmap["score"] = "BETA"
varmap["SNP"] = "SNP"


### add a menu ###

nvar=len(sys.argv)
cvar=1

while nvar - cvar > 2:
    tag=sys.argv[cvar]
    if tag == "--ea":
        cvar+=1
        varmap["EA"]=sys.argv[cvar]
    elif tag == "--pv":
        cvar+=1
        varmap["pvalue"]=sys.argv[cvar]
    elif tag == "--score":
        cvar+=1
        varmap["score"]=sys.argv[cvar]
    elif tag == "--snpname":
        cvar+=1
        varmap["SNP"]=sys.argv[cvar]
    else:
        sys.stderr.write("ERROR. Unrecognized option:" + tag + "\n")
        help(sys.argv[0], varmap)
    cvar += 1

if nvar - cvar < 2:
    help(sys.argv[0], varmap)

inname = sys.argv[cvar]
outbase = sys.argv[cvar+1]

processSumStat(inname, outbase, varmap)
