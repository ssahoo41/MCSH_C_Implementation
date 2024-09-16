import os
import pickle
import lmdb
import numpy as np
import ase.io
import torch
from tqdm import tqdm
from amptorch.preprocessing import AtomsToData, FeatureScaler, TargetScaler
from amptorch.descriptor.GMPOrderNorm import GMPOrderNorm
from ase import Atoms
from ase.calculators.emt import EMT
from json import load
from ase.calculators.singlepoint import SinglePointCalculator
import sys

from amptorch.ase_utils import AMPtorch
from amptorch.trainer import AtomsTrainer


def construct_lmdb(images, lmdb_path="./data.lmdb", normaliers_path="./normalizers.pt",order=0):
    """
    images: list of ase atoms objects (or trajectory) for fingerprint calculatation
    lmdb_path: Path to store LMDB dataset.
    normaliers_path: path of the scalers, create and store them if not exist
    """
    db = lmdb.open(
        lmdb_path,
        map_size=1099511627776 * 2,
        subdir=False,
        meminit=False,
        map_async=True,
    )

    # Define GMP
    sigmas = np.logspace(-1,0.6,3).round(2)


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
    MCSHs = MCSHs_dict[order]

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

    fp_params = GMP

    training_atoms = images
    elements = np.array([atom.symbol for atoms in training_atoms for atom in atoms])
    elements = np.unique(elements)
    descriptor = GMPOrderNorm(MCSHs=fp_params,elements=elements)
    descriptor_setup = ("gmpordernorm", GMP, GMP.get('cutoff'), elements)
    forcetraining = False

    a2d = AtomsToData(
        descriptor=descriptor,
        r_energy=True,
        r_forces=False,
        save_fps=False,
        fprimes=forcetraining,
    )

    data_list = []
    idx = 0
    for image in tqdm(
        images, desc="calculating fps", total=len(images), unit=" images",
    ):
        do = a2d.convert(image, idx=idx)
        data_list.append(do)
        idx += 1

    if os.path.isfile(normaliers_path):
        normalizers = torch.load(normaliers_path)
        feature_scaler = normalizers["feature"]
        target_scaler = normalizers["target"]

    else:
        scaling = {"type": "normalize", "range": (-1, 1),"elementwise":False}
        feature_scaler = FeatureScaler(data_list, forcetraining, scaling)
        target_scaler = TargetScaler(data_list, forcetraining)
        normalizers = {
            "target": target_scaler,
            "feature": feature_scaler,
        }
        torch.save(normalizers, normaliers_path)

    feature_scaler.norm(data_list)
    target_scaler.norm(data_list)

    idx = 0
    for do in tqdm(data_list, desc="Writing images to LMDB"):
        txn = db.begin(write=True)
        txn.put(f"{idx}".encode("ascii"), pickle.dumps(do, protocol=-1))
        txn.commit()
        idx += 1

    txn = db.begin(write=True)
    txn.put("feature_scaler".encode("ascii"), pickle.dumps(feature_scaler, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put("target_scaler".encode("ascii"), pickle.dumps(target_scaler, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put("length".encode("ascii"), pickle.dumps(idx, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put("elements".encode("ascii"), pickle.dumps(elements, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put(
        "descriptor_setup".encode("ascii"), pickle.dumps(descriptor_setup, protocol=-1)
    )
    txn.commit()

    db.sync()
    db.close()


if __name__ == "__main__":
    torch.set_default_tensor_type(torch.DoubleTensor)
    images = []

    dataset_name = 'troubleshoot_dataset.json'
    order = 4
    folder_prefix = 'test_gmpon'
    
    traj_train, e_train, CHO_count_train, traj_test, e_test, CHO_count_test  = load(open('/home/lucas/pace_hive/data/NNFF_calc/paper1_data/DFT/'+dataset_name))#should use your own path here
    store_folder = './'+folder_prefix
    if not os.path.exists(store_folder):
        os.mkdir(store_folder)
    os.chdir(store_folder)

    images = []
    for traj,e in zip(traj_train, e_train):
        fn = ase.io.read('/home/lucas/local_data/ggusmao/DFT/'+traj)
        force = fn.get_forces()
        calc = SinglePointCalculator(fn, energy=e, forces=force)
        fn.set_calculator(calc)
    
        images.append(fn)

    print('load train!')

    construct_lmdb(images, lmdb_path="./data_testgmpon.lmdb", normaliers_path="./normalizers_testglmpon.pt",order=order)