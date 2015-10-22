#import nose
import os
import sys
import numpy as np
sys.path.append(".")

from lintools import extract_data, visualisation

a = extract_data.ContactsData()

a.from_simdata("data/docked_ami_100ns_nowater.gro", "data/docked_ami_50ns_1_100_whole_no_water_skip100.xtc","protein","not protein")

def test_imported_shape_correct():
	assert a.distancematrix.shape == (101, 19564, 45)

def test_imported_values_correct():
	assert (a.distancematrix < 5).sum() == 69795

a.to_cache("deleteme_py_test_data")

def test_saved_all_files():
	assert os.path.exists("deleteme_py_test_data/distances.npy")
	assert os.path.exists("deleteme_py_test_data/ligand.json")
	assert os.path.exists("deleteme_py_test_data/protein.json")

b = extract_data.ContactsData()
b.from_cache("deleteme_py_test_data")

def test_loaded_shape_correct():
	assert b.distancematrix.shape == (101, 19564, 45)

def test_loaded_values_correct():
	assert (b.distancematrix < 5).sum() == 69795

def test_resids_correct():
	assert list(b.protein_resids())[::200] == [1, 12, 24, 38, 51, 65, 78, 90, 104, 116, 127, 140, 152, 166, 179, 190, 204, 216, 231, 243, 255, 267, 280, 292, 306, 320, 335, 346, 359, 371, 383, 396, 411, 423, 436, 448, 460, 473, 485, 497, 511, 522, 537, 550, 562, 576, 588, 600, 615, 629, 640, 652, 665, 676, 689, 702, 714, 727, 738, 752, 762, 773, 787, 802, 815, 827, 838, 852, 865, 878, 890, 902, 914, 927, 939, 952, 967, 979, 991, 1006, 1018, 1031, 1046, 1057, 1070, 1081, 1093, 1108, 1120, 1134, 1146, 1159, 1171, 1185, 1198, 1210, 1222, 1235]

def test_nearest_atoms_correct():
	rdict = b.nearest_atoms(310)[0]
	assert np.abs(rdict['std'] - 1.3396699) < 0.0001
	assert rdict['name'] == '1683,H06'
	assert np.abs(rdict['mean'] - 3.8543625) < 0.0001




