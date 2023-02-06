import os
import h5py as hpy
import numpy as np
from tqdm import tqdm

naked=np.zeros((1,15))
#path = '/home/lucas/githubRepos/MCSH_C_Implementation/GMP/processed/descriptors/GMPOrderNorm/64a9acbbcc41f61d4cf501fbe23eadd1'
path = '/home/lucas/githubRepos/MCSH_C_Implementation/processed/descriptors/GMPOrderNorm/64a9acbbcc41f61d4cf501fbe23eadd1'
os.chdir(path)
descriptor_files = os.listdir()
for file in tqdm(descriptor_files,total=len(descriptor_files)):
    file_path = os.path.join(path,file)
    try:
        f = hpy.File(file_path,'r')
    except:
        print('{} not loaded properly.\n'.format(file))
        continue
    dset = f['0']
    elem_pres = dset.keys()
    for elem in elem_pres:
        j = dset[elem]['fps'].shape[0]
        for k in range(j):
            naked = np.vstack((naked,dset[elem]['fps'][k,:]))
'''           
for elem in vals_dict.keys():
    vals_dict[elem] = np.concatenate(vals_dict[elem],axis=0)
    final_holder[elem].append(vals_dict[elem])
'''
print('Loading ccydata completed completed')
naked = naked[1:,:]
print(np.mean(naked,axis=0))
print(np.std(naked,axis=0))
