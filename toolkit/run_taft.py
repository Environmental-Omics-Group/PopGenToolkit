#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
TAFT Parallel Runner

author(s): Albert (aki) Zhou
added: 14-09-2020

"""


import os
import numpy as np
from kagami.comm import available, smap, unpack, filePath, filePrefix, fileTitle, checkInputFile, checkOutputFile
from kagami.portals import tablePortal
from kagami.wrappers import BinaryWrapper


def _slb(index, batch, nbatches, binpath, nprocs):
    print(f'  > loading index from [{index}]')
    ifiles = np.array(tablePortal.loadcsv(index)).flatten()
    print(f'  > [{len(ifiles)}] file listed')

    if available(batch): 
        ifiles = np.array_split(ifiles, nbatches)[batch-1]
        print(f'  > running batch [{batch} of {nbatches}] with [{len(ifiles)}] inputs')
    else:
        print(f'  > running all batches with [{len(ifiles)}] inputs')
    print(f'  > info files [{ifiles[0]} - {ifiles[-1]}]')

    _ipath = lambda x: os.path.join(filePath(index), x)
    params = [(_ipath(fn), _ipath(fileTitle(fn) + '_out.csv')) for fn in ifiles]
    smap(params, unpack(lambda ifn,ofn: (checkInputFile(ifn), checkOutputFile(ofn))))
    print(f'  > file paths verified')

    print(f'  > now running ...')
    bwp = BinaryWrapper(binpath)
    res = bwp.mapexec(params, nprocs = nprocs)
    for (_, ofn), (rcode,rstrs) in zip(params, res):
        tablePortal.savecsv(np.array(rstrs + [f'return code = {rcode}']).reshape((-1,1)), filePrefix(ofn) + '.log')

    print(f'  > done')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('index', help = 'path of the index file')
    parser.add_argument('-b', '--batchid', help = 'number of the current batch', type = int, default = None)
    parser.add_argument('-n', '--batches', help = 'total number of batches', type = int, default = 10)
    parser.add_argument('-x', '--exepath', help = 'TAFT executable path', default = './bin/TAFT-G')
    parser.add_argument('-p', '--procs', help = 'number of processors for each batch', type = int, default = 20)
    args = parser.parse_args()

    _slb(args.index, batch = args.batchid, nbatches = args.batches, binpath = args.exepath, nprocs = args.procs)
