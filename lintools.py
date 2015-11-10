if __name__ == '__main__':
	from ligand_interaction_network import LigandInteractions 
	from argparse import ArgumentParser
	import os

	#################################################################################################################

	parser = ArgumentParser(description='Counts the occurance of protein residues within a specified cutoff from a ligand from a simulation file.')
	#parser.add_argument('-j', '--jobname', dest = 'jobname', help ='Name of folder in which output files will be stored')
	parser.add_argument('-i', '--topology', dest = 'grofile', help='Input File name of topology file. Accepts gro, pdb files')
	parser.add_argument('-t', '--trajectory', dest = "xtcfile", nargs="*", default='None', help='Input File name of trajectory file(s). Accepts up to 3 xtc files (Optional. Default: "None")')
	parser.add_argument('-o', '--outname', dest = "output_name", help='Name of the output svg file.')
	#parser.add_argument('-p', '--input4', dest = "prot_sel", default = 'protein', help='MDAnalysis string for selection of protein (Optional. Default: "protein").')
	parser.add_argument('-c', '--cutoff', dest = "cutoff", default = 3.5, help='Input cutoff distance from the ligand that is taken into account in angstroms (Example: 3.5).')
	parser.add_argument('-ro', '--residueoffset', dest = "offset", default = 0, help='Input the number of offset residues for the protein. (Optional, default is 0)')
	parser.add_argument('-d', '--diagram_type', dest = "diagram_type", default="amino", help='Input type of diagramm required. Options: "clock" for clock diagrams (Only available with trajectory present), "helices" (NOT YET IMPLEMENTED) for diagrams representing residue membership to certain helices (requires user input to determine helices), "amino" showing the amino acid type')
	parser.add_argument('-hf', '--helix_file', dest = "helix_file", default="None", help='Input file for helices of your protein. To see the required format, check README or our GitHub page')
	args = parser.parse_args()

	###################################################################################################################
	assert len(args.output_name)>0, "Provide a name for output file"
	assert len(args.grofile)>0, "Provide input file"
	import MDAnalysis
	gro = MDAnalysis.Universe(args.grofile)
	list_of_non_ligands=["HOH","ARG","LYS","HIS","ASP","GLU","SER","THR", "ASN","GLN","PHE","TYR","TRP","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET"]
	potential_ligands={}
	i=0
	for resname in gro.residues.resnames:
	    if resname not in list_of_non_ligands and resname not in potential_ligands.values():
	        potential_ligands[i]=resname
	        i+=1
	for lig in potential_ligands:
	    print lig, potential_ligands[lig]
	ligand_name=potential_ligands[int(raw_input( "Choose a ligand to analyse:"))]

	lig = LigandInteractions()
	lig.load_universe(args.grofile, args.offset)
	assert args.diagram_type == "clock"  or args.diagram_type== "helices" or args.diagram_type=="amino", "Invalid diagramm type, try amino, helices or clock"
	if args.xtcfile=='None':
		assert args.diagram_type!="clock", "It is not possible to analyse data without trajectory, change your diagram type" 
		#can make diagramms of type "amino" or "helices"
	else:
		lig.occurance_count(args.grofile, args.xtcfile, ligand_name, args.cutoff, args.offset)
	if args.diagram_type=="helices":
		assert args.helix_file!="None", "A helix input file is required to plot a diagramm with helices"
		
	lig.get_closest_residues(ligand_name, args.cutoff)
	lig.get_closest_ligand_atoms(ligand_name)
	if args.diagram_type=="clock":
		lig.plot_clock_diagramms()
	elif args.diagram_type=="helices":
		assert args.helix_file!="None", "A helix input file is required to plot a diagramm with helices"
		lig.define_helices(args.helix_file, args.offset)
		lig.plot_helices_diagramms()
	elif  args.diagram_type=="amino":
		lig.define_amino_acids()
		lig.plot_amino_diagramms()

	lig.work_on_ligand(ligand_name)
	lig.convex_hull()
	lig.make_new_projection_values()
	lig.get_projection_values()
	lig.get_x_y_values()
	lig.draw_figure(args.diagram_type, ligand_name, args.output_name, args.grofile, args.xtcfile)
	file_list=[ligand_name+".pdb", "molecule_modified1.svg", "molecule_final.svg", "molecule.svg"]
	for f in file_list:
		os.remove(f)
	print "Ready!"

	print "The figure you created has been saved under filename "+args.output_name+".svg"
	print "This file format can be converted to any other image file format using "
	print "free software such as ImageMagick, Inkscape or GIMP."
	print "    "
	print "#########################################################################"
	print "       "
	print "Thank you for using Ligand Interaction Network. Check our GitHub account"
	print "(insert link). Developers: Laura Domicevica, Tom Newport and Teresa Paramo"
	print "at SBCB, Oxford University."
