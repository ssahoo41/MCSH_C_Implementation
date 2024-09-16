import time
import numpy as np
import os
import lmdb
import pickle
import torch
import random
import sys
import torch

start = time.time()
home_dir = "/home/lucas/githubRepos/MCSH_C_Implementation"

lmdb_path = home_dir+"/GMP/test_gmpon/"+"data_testgmpon.lmdb"
full_path = lmdb_path

db_env = lmdb.open(
    full_path,
    subdir=False,
    readonly=True,
    lock=False,
    readahead=False,
    meminit=False,
    max_readers=1,
)
count = 0
with db_env.begin(write=False) as txn:
    for key, value in txn.cursor():
        _object = pickle.loads(txn.get(key))
        if count == 0:
            _fprints = _object.fingerprint.detach().numpy()
            _elems = _object.atomic_numbers.detach().numpy()
            count += 1
        else:
            try:
                _fprints = np.concatenate((_fprints,_object.fingerprint.detach().numpy()),axis=0)
                _elems = np.concatenate((_elems,_object.atomic_numbers.detach().numpy()),axis=0)
            except AttributeError:
                pass

#_indices = list(np.arange(_fprints.shape[0]))
#num_2_keep = int(0.2*_fprints.shape[0])
#_keep = random.sample(_indices,num_2_keep)
fprints = _fprints#[_keep,:]
elems = _elems#[_keep]
np.savetxt(home_dir+"/GMP/test_fprints.csv",fprints,delimiter=",")
np.savetxt(home_dir+"/GMP/test_elems.csv",elems,delimiter=",")