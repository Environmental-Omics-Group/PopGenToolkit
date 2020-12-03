#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
Generate TAFT Batch Inputs

author(s): Albert (aki) Zhou
added: 05-11-2020

"""


import os
import numpy as np
from operator import itemgetter
from kagami.comm import smap, pmap, fold, unpack, paste, checkall, checkany, fileName, checkInputFile, checkOutputDir
from kagami.portals import tablePortal
from kagami.dtypes import Table


def _loadCounts(fname):
    dm = np.array(tablePortal.loadtsv(fname)[1:])
    scf, pos, cnt0, cnt1 = dm.T[[0,1,-2,-1]]
    otab = Table(
        np.vstack([cnt0, cnt1]).T.astype(int),
        rownames = smap(zip(scf,pos), lambda x: paste(x, sep = '_')),
        colnames = ['count_0', 'count_1'],
        rowindex = {'scaffold': scf, 'pos': pos.astype(int)}
    )
    return otab


def _genInfo(tabs):
    sids = np.array(tabs[0].rows_)
    assert checkall(tabs[1:], lambda x: np.all(x.rows_ == sids))

    olns = smap(
        zip(sids, *smap(tabs, lambda x: x.X_)), 
        lambda x: [[x[0], 2]] + \
                  smap(enumerate(x[1:]), unpack(lambda i,v: [f'Alelle_freq_Sample_{i}'] + list(v))) + \
                  [[]] # empty line
    )
    return olns


def _genBatch(outlns, nbatch):
    hdptn = '''
Simulation_Test_choice(T|F)   T
Algorithm(1=EmpBayes|2=Full)  2
WapTest_choice(T|F)           T
ChiTest_choice(T|F)           T
Ploidia                       2
Sampling_plan(1|2)            1

Nr_of_temporal_samples        2
Nr_gen_sample_1-2             6

Variable_Ne(T|F)              F
N_value                       1660001
Ne_value                      1660000

Nr_of_simulations             2000000
Random_number_seed            -1

Nr_of_loci                    %d

Maximum_nr_of_alleles         2
    '''    
    _hds = lambda nloc: smap((hdptn % nloc).lstrip().split('\n'), lambda x: x.split())
    
    oids = np.array_split(np.arange(len(outlns)), nbatch)

    def _merge(ols):
        odm = ols[0]
        for d in ols[1:]: odm.extend(d)
        return odm
    odlst = smap(oids, lambda x: itemgetter(*x)(outlns), lambda x: _hds(len(x)) + _merge(x))

    return odlst
 

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('counts', help = 'vcf allele counts', nargs = '+')
    parser.add_argument('outpath', help = 'output folder path')
    parser.add_argument('-n', '--batches', help = 'total number of batches', type = int, default = 1000)
    args = parser.parse_args()

    # main
    ifiles = args.counts
    smap(ifiles, checkInputFile)
    print(f'  > [{len(ifiles)}] count files listed')

    ctabs = pmap(ifiles, _loadCounts)
    for f,t in zip(ifiles,ctabs): print(f'  > [{t.nrow}] loci loaded from [{fileName(f)}]')
        
    sids = fold(smap(ctabs, lambda x: x.rows_), lambda x,y: np.intersect1d(x,y))
    print(f'  > [{len(sids)}] loci shared in input files')

    ctabs = smap(ctabs, lambda x: x[sids])
    if checkany(ctabs, lambda x: 0 in x.X_):
        print('  > zero count(s) found, add 1 to all counts for valid TAFT input')
        for t in ctabs: t.X_ += 1
    
    print(f'  > generating batches ...')
    olns = _genInfo(ctabs)
    ochk = _genBatch(olns, args.batches)
    
    opath = args.outpath
    checkOutputDir(opath)
    print(f'  > saving info files to [{opath}]')
        
    ofiles = np.array([f'info_chunk{i}.dat' for i in range(len(ochk))])
    for f,o in zip(ofiles,ochk): tablePortal.save(o, fname = os.path.join(opath,f), delimiter = ' ')

    print(f'  > saving index file')
    tablePortal.savecsv(ofiles.reshape((-1,1)), os.path.join(opath, 'index.csv'))

    print(f'  done.')
