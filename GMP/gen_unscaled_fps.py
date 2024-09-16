from amptorch.descriptor.GMPOrderNorm import GMPOrderNorm
from amptorch.preprocessing import (
    AtomsToData,
    FeatureScaler,
    TargetScaler
)
import pickle
import sys
from json import load
from ase.io import read
import os
from ase.calculators.singlepoint import SinglePointCalculator
import torch
import numpy as np

#Implement ability to use established scaling for test set - load feature scaler from normalizer.pt
def process(a2d,images,forcetraining,scaling,testset=False):
    data_list = a2d.convert_all(images)
    if testset:
        normalizers = torch.load('./normalizers.pt')
        feature_scaler = normalizers['feature']
        target_scaler = normalizers['target']
        data_list = feature_scaler.norm(data_list)
        data_list = target_scaler.norm(data_list)
    else:
        feature_scaler = FeatureScaler(data_list,forcetraining,scaling)
        target_scaler = TargetScaler(data_list,forcetraining)
        data_list = feature_scaler.norm(data_list)
        data_list = target_scaler.norm(data_list)

        normalizers = {'target': target_scaler,'feature':feature_scaler}
        torch.save(normalizers,'normalizers.pt')

    return data_list

sigmas = sigmas = np.logspace(-1,0.6,3).round(2)

MCSHs_dict = {
0: { "orders": [0], "sigmas": sigmas,},
1: { "orders": [0,1], "sigmas": sigmas,},
2: { "orders": [0,1,2], "sigmas": sigmas,},
3: { "orders": [0,1,2,3], "sigmas": sigmas,},
4: { "orders": [0,1,2,3,4], "sigmas": sigmas,},
5: { "orders": [0,1,2,3,4,5], "sigmas": sigmas,},
6: { "orders": [0,1,2,3,4,5,6], "sigmas": sigmas,},
7: { "orders": [0,1,2,3,4,5,6,7], "sigmas": sigmas,},
8: { "orders": [0,1,2,3,4,5,6,7,8], "sigmas": sigmas,},
9: { "orders": [0,1,2,3,4,5,6,7,8,9], "sigmas": sigmas,},
}
#The atom_gaussians should be modified by where you save your own Gaussians, if you do not have those Gaussian potentials, please check /storage/hive/project/chbe-medford/medford-share/shared/cchang373/valence_gaussians
MCSHs = MCSHs_dict[4]
'''
GMP = { "MCSHs": MCSHs,
        "atom_gaussians": {
            "C": "/home/lucas/local_data/pseudo/C_pseudodensity_4.g",
            "O": "/home/lucas/local_data/pseudo/O_pseudodensity_4.g",
            "H": "/home/lucas/local_data/pseudo/H_pseudodensity_2.g",
            "Pt": "/home/lucas/local_data/pseudo/Pt_pseudodensity_6.g",
            "Au": "/home/lucas/local_data/pseudo/Au_pseudodensity_6.g",
            "Rh": "/home/lucas/local_data/pseudo/Rh_pseudodensity_4.g",
            "Pd": "/home/lucas/local_data/pseudo/Pd_pseudodensity_4.g",
            "Ag": "/home/lucas/local_data/pseudo/Ag_pseudodensity_4.g",
        },
        "cutoff": 18,
        "solid_harmonics":False
    }
'''
GMP = { "MCSHs": MCSHs,
        "atom_gaussians": {
            "Al": "/home/lucas/local_data/pseudo/Al_pseudodensity_4.g",
        },
        "cutoff": 18,
        "solid_harmonics":False
    }

fp_params = GMP
forcetraining = True

#only need this if atoms.get_potential_energy wouldn't work - need to attach energy
#Probably depends on data set - assume I need to for ccy data - could be eliminated with preprocessing step
'''
dataset_name = 'troubleshoot_dataset.json'


traj_train, e_train, CHO_count_train, traj_test, e_test, CHO_count_test = load(open('/home/lucas/pace_hive/data/NNFF_calc/paper1_data/DFT/'+dataset_name))#should use your own path here
'''
traj_train, e_train = ['Al_scell.traj'],[0.]
folder_prefix = 'test_gmpon'
store_folder = './'+folder_prefix
if not os.path.exists(store_folder):
    os.mkdir(store_folder)
os.chdir(store_folder)
#train and test images
images_train, images_test = [], []
elements = ['Al']

for traj,e in zip(traj_train, e_train):
    fn = read('/home/lucas/githubRepos/MCSH_C_Implementation/'+traj)
    force = 0.
    #force = fn.get_forces()
    calc = SinglePointCalculator(fn, energy=e, forces=force)
    fn.set_calculator(calc)
    
    images_train.append(fn)

print('load train!')
'''
for traj,e in zip(traj_test, e_test):
    fn = read('/storage/hive/project/chbe-medford/medford-share/shared/cchang373/DFT/'+traj)
    force = fn.get_forces()
    calc = SinglePointCalculator(fn, energy=e, forces=force)
    fn.set_calculator(calc)

    images_test.append(fn)

print('load test!')
'''
scaling = {'type':'normalize','range':(-1,1),'elementwise':False}

descriptor = GMPOrderNorm(MCSHs=fp_params,elements=elements)

a2d = AtomsToData(
    descriptor=descriptor,
    r_energy=True,
    r_forces=forcetraining,
    save_fps=True,
    fprimes=forcetraining,
    cores=1
)
testset=False
train_list = process(a2d,images_train,forcetraining,scaling,testset=testset)

save_obj = False
if save_obj:
    if not testset:
        with open('{}_scaled_desc'.format(folder_prefix),'wb') as fp:
            pickle.dump(train_list,fp)
    else:
        with open('{}_scaled_desc_test'.format(folder_prefix),'wb') as fp:
            pickle.dump(train_list,fp)