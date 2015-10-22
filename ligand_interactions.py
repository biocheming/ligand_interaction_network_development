""" Script for extracting distance/contact data from simulation files and directly visualising.
Raw distance data is saved as {jobname}/distances.npy 
        :Arguments:
            *grofile*
                string of the name of the GROMACS gro file
            *xtcfile*
                trajectory file
            *ligand_name*
                string of the ligand name
            *cutoff*
                cutoff distance from the ligand that is taken into account in angstroms
        :Returns:

"""

import MDAnalysis as mda
import numpy as np
from sklearn.metrics import pairwise
import os
import json
import matplotlib.pyplot as plt
import sys
from argparse import ArgumentParser

import Visualisation
import ExtractData

##########################################################################

parser = ArgumentParser(description='Counts the occurance of protein residues within a specified cutoff from a ligand from a simulation file.')
parser.add_argument('-j', '--jobname', dest = "jobname", help ='Name of folder in which output files will be stored')
parser.add_argument('-i', '--input1', dest = "topfile", help='Input File name of topology file (Example: docked_ami_100ns_nowater.gro).')
parser.add_argument('-t', '--input2', dest = "trajfile", help='Input File name of trajectory file (Example: docked_ami_50ns_1_100_whole_no_water_skip100.xtc).')
parser.add_argument('-l', '--input3', dest = "lig_sel", default='not protein', help='MDAnalysis string for selection of ligand (Optiona. Default: "not protein").')
parser.add_argument('-p', '--input4', dest = "prot_sel", default = 'protein', help='MDAnalysis string for selection of protein (Optional. Default: "protein").')
parser.add_argument('-d', '--input5', dest = "distcutoff", default = 5, help='NOT YET IMPLEMENTED. Input cutoff distance from the ligand that is taken into account in angstroms (Example: 5).')
args = parser.parse_args()

##########################################################################

# Measure (and save) the distances
found_contacts = ExtractData.ContactsData()
found_contacts.from_simdata(args.topfile, args.trajfile, args.prot_sel, args.lig_sel)
found_contacts.to_cache(args.jobname)

# Do some analysis
x = sorted(list(set(found_contacts.protein_resids())))
y = found_contacts.ligand_atoms
Visualisation.sort_of_works(x, y, found_contacts)



