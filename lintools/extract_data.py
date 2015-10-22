""" Module containing the ContactsData class for measuring protein residue - ligand atom distances,
saving and reading these distances, and some level of manipulation of the data.

If run as a script, the contacts in a supplied gro/trajectory file are measured and saved. 
"""

import MDAnalysis as mda
import numpy as np
from sklearn.metrics import pairwise
import os
import json
import matplotlib.pyplot as plt
import sys

##########################################################################

class ContactsData:
	"""For the storing and reading of contacts between a protein and ligand"""
	def __init__(self):
		self.protein_atoms = []
		self.ligand_atoms = []
		self.distancematrix = None
	def from_simdata(self, structure, trajectory, protein_selector, ligand_selector):
		""" Extract contact and atom name data from simulation and populate class attribute """
		universe = mda.Universe(structure, trajectory)
		protein = universe.select_atoms(protein_selector)
		ligand = universe.select_atoms(ligand_selector)
		frames = []
		for ts in universe.trajectory:
			frames.append(pairwise.pairwise_distances(protein.atoms.positions, ligand.atoms.positions))
		frames = np.array(frames)
		self.protein_atoms = [str(x.resid) + "," + x.name for x in protein.atoms]
		self.ligand_atoms = [str(x.resid) + "," + x.name for x in ligand.atoms]
		# Frame, Protein Atom, Ligand Atom
		self.distancematrix = frames
	def to_cache(self, file_location):
		""" Save class attributes to disk"""
		if not os.path.exists(file_location):
			os.makedirs(file_location)
		np.save(os.path.join(file_location, "distances.npy"), self.distancematrix)
		with open(os.path.join(file_location,"protein.json"),"w") as fh:
			json.dump(self.protein_atoms, fh)
		with open(os.path.join(file_location,"ligand.json"),"w") as fh:
			json.dump(self.ligand_atoms, fh)
	def from_cache(self, file_location):
		""" Populate class attributes from disk """
		for in_file in ("distances.npy", "protein.json", "ligand.json"):
			if os.path.exists(os.path.join(file_location, in_file)):
				pass
			else:
				raise IOError("File {name} not found!".format(name = in_file))
		self.distancematrix = np.load(os.path.join(file_location, "distances.npy"))
		with open(os.path.join(file_location,"protein.json"),"r") as fh:
			self.protein_atoms = json.load(fh)
		with open(os.path.join(file_location,"ligand.json"),"r") as fh:
			self.ligand_atoms = json.load(fh)
	def protein_resids(self):
		""" Returns list of protein residue ids (for each atom) """
		return np.array([int(x.split(",")[0]) for x in self.protein_atoms])
	def filter_by_residue(self, residue_id):
		""" Returns contact data for atoms in selected residue """
		protein_selector = self.protein_resids() == residue_id
		# Frame, Ligand Atom
		return self.distancematrix[:,protein_selector,:]
	def contact_time(self, residue_id, distance):
		""" Returns array of existence of a contact within the cuttoff distance for the selected residue """
		return (self.filter_by_residue(residue_id)<distance).any(axis=1).any(axis=1)
	def contact_times(self, cutoff = 5):
		""" Return dictionary of total number of frames spent in contact with ligand for each residue """
		return {resid:self.contact_time(resid, cuttoff).sum() for resid in list(set(self.protein_resids()))}
	def protein_distances(self, residue_id):
		""" Returns array of minimum distance of selected residue to each ligand atom in each frame """
		return self.filter_by_residue(residue_id).min(axis=1)
	def nearest_atoms(self, residue_id, cutoff=5):
		""" For the selected residue, returns ranked list of average distance to each ligand atom """
		distance_data = self.protein_distances(residue_id)
		means = distance_data.mean(axis=0)
		stds = distance_data.std(axis=0)
		relevant_atoms = means < cutoff
		r_means = means[relevant_atoms]
		r_stds = stds[relevant_atoms]
		r_data = np.array(self.ligand_atoms)[relevant_atoms]
		return sorted([{"name":d, "std":s, "mean":m} for m, s, d in zip(r_means, r_stds, r_data)], key= lambda x : x["mean"])

##########################################################################
		
if __name__ == "__main__":
	# Set up options for when run as script...
	from argparse import ArgumentParser
	parser = ArgumentParser(description='Calculates all protein-ligand distances and saves to {jobname}/distances.npy')
	parser.add_argument('-j', '--jobname', dest = "jobname", help ='Name of folder in which output files will be stored')
	parser.add_argument('-i', '--input1', dest = "topfile", help='Input File name of topology file (Example: docked_ami_100ns_nowater.gro).')
	parser.add_argument('-t', '--input2', dest = "trajfile", help='Input File name of trajectory file (Example: docked_ami_50ns_1_100_whole_no_water_skip100.xtc).')
	parser.add_argument('-l', '--input3', dest = "lig_sel", default='not protein', help='MDAnalysis string for selection of ligand (Optional. Default: "not protein").')
	parser.add_argument('-p', '--input4', dest = "prot_sel", default = 'protein', help='MDAnalysis string for selection of protein (Optional. Default: "protein").')
	args = parser.parse_args()
	
	# Read all contacts from the provided files and save out as appropriate
	found_contacts = ContactsData()
	found_contacts.from_simdata(args.topfile, args.trajfile, args.prot_sel, args.lig_sel)
	found_contacts.to_cache(args.jobname)
	
