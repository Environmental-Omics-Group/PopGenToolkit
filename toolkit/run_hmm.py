#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
TAFT Parallel Runner

author(s): Albert (aki) Zhou
added: 11-10-2019

"""


import os
import numpy as np
from kagami.comm import available, smap, filePath, fileTitle, checkInputFile
from kagami.portals import tablePortal
from kagami.wrappers import BinaryWrapper


def _slb(index, batch, nbatches, rpath, script, nprocs):
    print(f'  > loading file list from [{index}]')
    ifiles = np.array(tablePortal.loadcsv(index)).flatten()
    print(f'  > [{len(ifiles)}] file listed')

    if available(batch): 
        ifiles = np.array_split(ifiles, nbatches)[batch-1]
        print(f'  > running batch [{batch} of {nbatches}] with [{len(ifiles)}] inputs')
    else:
        print(f'  > running all batches with [{len(ifiles)}] inputs')
    print(f'  > info files [{ifiles[0]} - {ifiles[-1]}]')
    
    _ipath = lambda x: os.path.join(filePath(index), x)
    params = [(script, _ipath(fn), nprocs) for fn in ifiles]
    smap(params, lambda x: checkInputFile(x[1]))
    print(f'  > file paths verified')

    print(f'  > now running ...')
    bwp = BinaryWrapper(rpath)
    res = bwp.mapexec(params)

    print('  > done')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('script', help = 'HMM script')
    parser.add_argument('index', help = 'index file path')
    parser.add_argument('-b', '--batchid', help = 'number of the current batch', type = int, default = None)
    parser.add_argument('-n', '--batches', help = 'total number of batches', type = int, default = 10)
    parser.add_argument('-r', '--rspath', help = 'RScript executable path', default = 'Rscript')
    parser.add_argument('-p', '--procs', help = 'number of processors for each batch', type = int, default = 20)
    args = parser.parse_args()

    _slb(args.index, args.batchid, args.batches, args.rspath, args.script, args.procs)

