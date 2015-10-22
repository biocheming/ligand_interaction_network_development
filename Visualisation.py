"""Module containing functions for various visualisations of protein-ligand distance data.

If run as a function, visualisation is performed using the saved contact data in the supplied directory

CURRENTLY ONLY CONTAINS A SIMPLE VISUALISATION OF MEAN DISTANCE FOR EACH PROTEIN/LIGAND.
"""

import MDAnalysis as mda
import numpy as np
from sklearn.metrics import pairwise
import os
import json
import matplotlib.pyplot as plt

import sys
import ExtractData

##########################################################################

def sort_of_works(x, y, contacts):
	""" Quick and dirty visualisation of mean residue/ligand-atom distance """
	data2D = np.zeros((len(x), len(y)))
	for residue in x:
		x_index = x.index(residue)
		atoms = contacts.nearest_atoms(residue, cutoff=100)
		for atom in atoms:
			y_index = y.index(atom["name"])
			data2D[x_index, y_index] = atom["mean"]
	plt.imshow(data2D, interpolation="nearest", cmap=plt.get_cmap("Spectral"))	
        plt.show()
        #plt.savefig('out.png')

#### MORE VISUALISATION STUFF...?

	
if __name__ == "__main__":
	# Set up options for when run as script...
	from argparse import ArgumentParser
	parser = ArgumentParser(description='Counts the occurance of protein residues within a specified cutoff from a ligand, given the residue-ligand distances in {name}/distances.npy.')
	parser.add_argument('-j', '--name', dest = "name", help ='Name of folder in which output files are stored')
	parser.add_argument('-d', '--input4', dest = "distcutoff", default = 5, help='NOT YET IMPLEMENTED. Input cutoff distance from the ligand that is taken into account in angstroms (Example: 5).')
	# options for various visualisations/output file names
	args = parser.parse_args()
	
	# Load the contact data
	contacts_loaded = ExtractData.ContactsData()
	contacts_loaded.from_cache(args.name)

	# Manipulate + visualise the data
	x = sorted(list(set(contacts_loaded.protein_resids())))
	y = contacts_loaded.ligand_atoms
	sort_of_works(x, y, contacts_loaded)

	### MORE VISUALISATION STUFF...?
