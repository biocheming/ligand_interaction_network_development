import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

from rdkit import Chem

from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from itertools import combinations

import matplotlib.cm as cmx
import matplotlib.colors as colors
import pylab

from collections import Counter
import json
import sys
import fileinput

from MDAnalysis import analysis
from MDAnalysis.analysis import distances

from shapely import geometry




class LigandInteractions(object):
    print "Testing for new residue dispersion"

    def renumber_residues(self, universe, offset):
        """Renumbers protein residues according to the offset value"""
        protein = universe.select_atoms("protein").residues
        protein.set_resids(protein.resids+int(offset))
        protein.resids


    def load_universe(self, grofile, offset):
        """Loads grofile as MDAnalysis universe"""
        self.topology = MDAnalysis.Universe(grofile)
        self.protein = self.topology.select_atoms("protein").residues
        self.protein.set_resids(self.protein.resids+int(offset))
        self.protein.resids
        self.occurance_run=False

    def occurance_count(self, grofile, xtcfile, ligand_name, cutoff, offset):
        """Counts the occurance of protein residues within a specified cutoff from a ligand and calculates which residues have been within cutoff for more than 50 ns (by default).
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
                *residue_counts*
                    a dictionary of residues and how many times during the simulation they are within
                    a cutoff value from the ligand"""
        i=0
        self.residue_counts={}
        for xtc in xtcfile:
            i+=1
            md_sim = MDAnalysis.Universe(grofile, xtc)
            lig.renumber_residues(md_sim, offset)
            frame_dict = {}
            firstframe_ps=None
           
            for frame in md_sim.trajectory:
                if ligand_name!="not protein":
                    selection = md_sim.select_atoms('protein and around '+str(cutoff)+' resname '+ligand_name)
                else:
                    selection = md_sim.select_atoms('protein and around '+str(cutoff)+' '+ligand_name)
                residue_list = [atom.resid for atom in selection]
                frame_dict[frame.time]=set(residue_list)

                if firstframe_ps == None:
                    firstframe_ps = frame.time
            
            #Calculate the cutoff time in ns - at the moment half of the simulation time
            lastframe_time = max([f for f in frame_dict.keys()])
            lastframe_time_ns = lastframe_time/1000
            self.cutoff_time_ns = int(lastframe_time_ns/2)

            self.residue_counts[i] = Counter([item for sublist in frame_dict.values() for item in sublist])
        
        self.occurance_run=True


    def get_closest_residues(self, ligand_name, cutoff):
        """Find the list of residues to be plotted using cutoff"""
        self.list_of_plotted_res=[]
        if self.occurance_run==True:
            #if only 1 traj is analysed
            if len(self.residue_counts)==1:
                for res in self.residue_counts[1].keys():
                    if self.residue_counts[1][res]>self.cutoff_time_ns:
                        self.list_of_plotted_res.append(res)
            else :
                new_res_list={}
                for xtc in self.residue_counts:
                    for res in self.residue_counts[xtc]:
                        new_res_list[res]=[]
                for res in new_res_list:
                    for xtc in self.residue_counts:
                        if res in self.residue_counts[xtc].keys():
                            new_res_list[res].append(self.residue_counts[xtc][res])
                        else:
                            new_res_list[res].append(0)
                for res in new_res_list:
                    for (index1, value1),(index2, value2) in combinations(enumerate(new_res_list[res]),2):
                        if value1>self.cutoff_time_ns and value2>self.cutoff_time_ns:
                            if res not in self.list_of_plotted_res:
                                self.list_of_plotted_res.append(res)
        else:
            if ligand_name!="not protein":
                selection = self.topology.select_atoms('protein and around '+str(cutoff)+' resname '+ligand_name)
            else:
                selection = self.topology.select_atoms('protein and around '+str(cutoff)+' '+ligand_name)
            residue_list = [atom.resid for atom in selection]
            for res in residue_list:
                if res not in self.list_of_plotted_res:
                    self.list_of_plotted_res.append(res)

        

    def get_closest_ligand_atoms(self, ligand_name):
        """Finds the ligand atom that is closest to a particular residue"""
        ## Selecting ligand without hydrogen atoms as these are not depicted in the RDKit fig
        if ligand_name!="not protein":
            self.ligand = self.topology.select_atoms("resname "+ligand_name+" and not name H*")
        else:
            self.ligand = self.topology.select_atoms(ligand_name+" and not name H*")
        lig_pos = self.ligand.positions
        self.closest_atoms = {}
        # Run the distance matrix calculation to get the closest atom (could be improved by using itertools combinations)
        for residue in self.list_of_plotted_res:
            residue_select= self.topology.select_atoms("resid "+str(residue))
            res_pos = residue_select.positions
            dist_array = MDAnalysis.analysis.distances.distance_array(lig_pos, res_pos)
            for atom in self.ligand:
                if dist_array[atom.id].min()== dist_array.min():
                    self.closest_atoms[residue]=atom.name, dist_array[atom.id].min()
                    #print "Residue "+str(residue)+" is closest to ligand atom "+atom.name+" with distance "+str(dist_array[atom.id].min())
        
        
        #Get the RES123 format
        renumbered_residues=[]
        for residue in range(len(self.protein.resids)):
            renumbered_residues.append(str(self.protein.resnames[residue])+str(self.protein.resids[residue]))
        #Get the newly formatted residues in a dictionary
        calculated_offset=self.protein.resids[0]
        self.dict_for_clock_diagramms={}
        if self.occurance_run==True:
            for res in self.list_of_plotted_res:
                self.dict_for_clock_diagramms[renumbered_residues[res-int(calculated_offset)]]=[]
            for xtc in self.residue_counts.values():
                for res in self.list_of_plotted_res:
                    self.dict_for_clock_diagramms[renumbered_residues[res-int(calculated_offset)]].append(xtc[res])
        else:
            for  res in self.list_of_plotted_res:
                self.dict_for_clock_diagramms[renumbered_residues[res-int(calculated_offset)]] = 1 
    #give an arbitrary value of 1 so that matplotlib has something to plot

        
    def define_helices(self, helices_text_file, offset):
        helices = {}
        self.colors_helices ={1:(0.5490,0.7961, 0.8),2:(0.7961, 0.3921568, 0.1843),3:(0.70588,0.36078, 0.82745),4:(0.52941, 0.831373, 0.262745),5:(0.337254, 0.21176, 0.37647),6:(0.34902, 0.45098, 0.20392),7:(0.32941, 0.21176, 0.145098),8:(0.796078, 0.647059, 0.58431),9:(0.53333, 0.545098, 0.8),10:(0.4941176, 0.82353, 0.545098),11:(0.784314, 0.32941,0.58431),12:(0.3098, 0.415686, 0.41961), 13:(0.6, 0.6, 0.6)}
            
        with open (helices_text_file, "r") as h:
            lines = h.readlines()
            for line in lines:
                helices[int(line.rsplit(",",3)[0])]=[int(line.rsplit(",",3)[1]),int(line.rsplit(",",3)[2])]

        i=0
        for res in helices:
            i+=1
            while i< len(helices):
                assert helices[res+1][0]-helices[res][1] > 0, "Helices are overlapping in your input file "+helices_text_file+". Please check!"
                break
    
        self.residues_within_helix={}
        for helix in helices:
            for res in self.dict_for_clock_diagramms:
                if int(res[3::])-int(offset)>=helices[helix][0] and int(res[3::])-int(offset)<=helices[helix][1]:
                    self.residues_within_helix[res]=helix
        print "Defining helices..."
    
    def define_amino_acids(self):
        #Group amino acids by charge and structure
        self.amino_acids = {"acidic":["ASP","GLU"], "basic":["LYS","ARG"], "aromatic":["PHE","TYR","TRP"],"polar":["SER","THR","ASN","GLN","CYS","HIS"],"hydrophobic":["ALA","VAL","ILE","LEU","MET","GLY"]}
        #Color scheme
        self.colors_amino_acids = {"acidic":"#889DCC", "basic":"#D06AC1", "aromatic":"#9FC74A", "polar":"#D9774B","hydrophobic":"#6AC297"}
        self.amino_acid_type={}
        for residue in self.dict_for_clock_diagramms.keys():
            for aa_type, aa in self.amino_acids.items():
                for amino_acid in aa:
                    if residue[0:3] == amino_acid:
                        self.amino_acid_type[residue]=aa_type
        print "Defining amino acid type..."

        #Matplotlib parameters
    matplotlib.rcParams['svg.fonttype'] = 'none'
    matplotlib.rcParams['font.family'] = ['sans-serif']
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['font.weight']=900
    matplotlib.rcParams['text.usetex'] = False
    matplotlib.rcParams['patch.linewidth'] = 0  

    def plot_clock_diagramms(self):
        """Uses matplotlib to plot clock diagramms used for data analysis of trajectories, for example, occurance time  over timecourse of simulation"""
        #Plot the residues in clock diagramm fashion
        colors_1=['#1f78b4','white']
        colors_2=['#33a02c','white']
        colors_3=['#6a3d9a','white']

        for res in self.dict_for_clock_diagramms.keys():
            plt.figure(figsize=(1.12, 1.12))
            if len(self.dict_for_clock_diagramms[res])==1:
                width=0.25
                ring1,_=plt.pie([self.dict_for_clock_diagramms[res][0],100-self.dict_for_clock_diagramms[res][0]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
            elif len(self.dict_for_clock_diagramms[res])==2:
                width=0.25
                ring1,_=plt.pie([self.dict_for_clock_diagramms[res][0],100-self.dict_for_clock_diagramms[res][0]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
                ring2,_=plt.pie([self.dict_for_clock_diagramms[res][1],100-self.dict_for_clock_diagramms[res][1]],  radius=1, startangle=90, colors=colors_2, counterclock=False)
            elif len(self.dict_for_clock_diagramms[res])==3:
                width=0.25
                ring1,_=plt.pie([self.dict_for_clock_diagramms[res][0],100-self.dict_for_clock_diagramms[res][0]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
                ring2,_=plt.pie([self.dict_for_clock_diagramms[res][1],100-self.dict_for_clock_diagramms[res][1]],  radius=1, startangle=90, colors=colors_2, counterclock=False)
                ring3,_=plt.pie([self.dict_for_clock_diagramms[res][2],100-self.dict_for_clock_diagramms[res][2]],  radius=1+width, startangle=90, colors=colors_3, counterclock=False)
            plt.axis('equal')
            if len(self.dict_for_clock_diagramms[res])==1:
                plt.setp(ring1, width=width)
                plt.text(0,-0.4,res[0:3]+"\n"+res[3::],ha='center',size=16, fontweight='bold')
            elif len(self.dict_for_clock_diagramms[res])==2:
                plt.setp(ring1+ring2, width=width)
                plt.text(0,-0.4,res[0:3]+"\n"+res[3::],ha='center',size=12, fontweight='bold')
            elif len(self.dict_for_clock_diagramms[res])==3:
                plt.setp(ring1+ring2+ring3, width=width)
                plt.text(0,-0.375,res[0:3]+"\n"+res[3::],ha='center',size=10, fontweight='bold')
            pylab.savefig(str(res[3::])+".svg", dpi=80, transparent=True)
            #plt.show()
        print "Plotting..."

    def plot_helices_diagramms(self):
        width=0.20        
        for res in self.dict_for_clock_diagramms:
            if res in self.residues_within_helix:
                color = [self.colors_helices[self.residues_within_helix[res]],'white']
                plt.figure(figsize=(1.12,1.12))
                ring1,_=plt.pie([1],  radius=1-width, startangle=90, colors=color, counterclock=False)
                plt.axis('equal')
                plt.setp(ring1, width=width, edgecolor='white')
                plt.text(0,-0.4,res[0:3]+"\n"+res[3::],ha='center',size=16, fontweight='bold')  
                pylab.savefig(str(res[3::])+".svg", dpi=80, transparent = True)
            else:
                color = [self.colors_helices[13], 'white'] # Matplotlib colors takes more than 1 color, so a second color must be added
                plt.figure(figsize=(1.12,1.12), dpi=100)
                ring1,_=plt.pie([1],  radius=1-width, startangle=90, colors=color, counterclock=False)
                plt.axis('equal')
                plt.setp(ring1, width=width, edgecolor='white')
                plt.text(0,-0.4,res[0:3]+"\n"+res[3::],ha='center',size=16, fontweight='bold')  
                pylab.savefig(str(res[3::])+".svg", dpi=80, transparent=True)
        print "Plotting..."

    def plot_amino_diagramms(self):
        for res in self.amino_acid_type:
            color = [self.colors_amino_acids[self.amino_acid_type[res]],'white']
            plt.figure(figsize=(1.12,1.12))
            #plot a random value (1) at the moment
            ring1,_=plt.pie([1],  radius=1, startangle=90, colors=color, counterclock=False)
            plt.axis('equal')
            plt.setp(ring1, width=1, edgecolor=self.colors_amino_acids[self.amino_acid_type[res]])
            plt.text(0,-0.5,res[0:3]+"\n"+res[3::],ha='center',size=18, fontweight="bold")
            pylab.savefig(str(res[3::])+".svg",dpi=80, transparent = True)
        
        print "Plotting..."


    def work_on_ligand(self, ligand_name):
        """Make initial SVG file of ligand that is later used to extract ligand atom coordinates"""
        # make ligand PDB file
        self.ligand.atoms.write(ligand_name+".pdb")
        # make MDA universe with just ligand
        ligand_mda = MDAnalysis.Universe(ligand_name+".pdb")
        #get atom names from MDanalysis and save in a list
        self.atom_names=[]
        for atom in ligand_mda.atoms.indices:
            self.atom_names.append(ligand_mda.atoms.names[atom])
        x=0
        while atom > len(self.atom_names)-1:
            assert ligand_mda.atoms.names[atom]!=ligand_mda.atoms.names[atom+1], 'Give your atoms unique names'
            x+=1

        # input the ligand in rkdit environment
        self.ligand_rdkit = Chem.MolFromPDBFile(ligand_name+".pdb")
        #Draw the initial SVG file from which the coordinates will be read
        AllChem.Compute2DCoords(self.ligand_rdkit)
        drawer = rdMolDraw2D.MolDraw2DSVG(600,300)
        opts = drawer.drawOptions()
        for i in range(self.ligand_rdkit.GetNumAtoms()):
            opts.atomLabels[i] = self.atom_names[i]
        drawer.DrawMolecule(self.ligand_rdkit)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:','')
        #write the svg figure to a file
        filesvg = open("molecule.svg", "w+")
        filesvg.write(svg)
        
    def convex_hull(self):
        """Draws a convex hull around ligand atoms and expands it, giving space to put diagramms on"""
        #Get coordinates of ligand atoms (needed to draw the convex hull around)
        ligand_atom_coords = []
        with open ("molecule.svg", "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("<text"): 
                    ligand_atom_coords.append([float(line.rsplit("'",6)[1])+100.00, float(line.rsplit("'",6)[3])+150.00]) #These numbers are added to make the ligand coordinates absolute rather than relative to their inner box (check SVG element)
        ligand_atom_coords=np.array(ligand_atom_coords)  
        # Get the convex hull around ligand atoms 
        self.a = geometry.MultiPoint(ligand_atom_coords).convex_hull

        # Get a dictionary with the coordinates of the *closest* atoms in the diagramm per residue
        self.x_coordinates_from_diagramm = {}
        self.y_coordinates_from_diagramm = {}
        with open ("molecule.svg", "r") as f:
            lines = f.readlines()
            for ligand_atom in self.closest_atoms:
                for line in lines:
                    if line.startswith("<text") and self.closest_atoms[ligand_atom][0] == line.rsplit("<",6)[2][6:] in line:
                        # add 100-150px to get absolute coordinates
                        self.x_coordinates_from_diagramm[ligand_atom] = float(line.rsplit("'",6)[1])+100.00
                        self.y_coordinates_from_diagramm[ligand_atom] = float(line.rsplit("'",6)[3])+150.00

        # Get the nearest point projection on the convex hull
        self.nearest_points_projection = {}
        nearest_points = {}
        self.nearest_points_coords ={}
        #Puts coords on different hulls
        for residue in self.x_coordinates_from_diagramm:
            b = self.a.boundary.parallel_offset(self.closest_atoms[residue][1]*32.0+24,"left",join_style=2).convex_hull
            point =geometry.Point((self.x_coordinates_from_diagramm[residue],self.y_coordinates_from_diagramm[residue]))
            self.nearest_points_projection[residue] = (b.boundary.project(point) % b.boundary.length)
            nearest_points[residue] = b.boundary.interpolate(self.nearest_points_projection[residue] % b.boundary.length)
            self.nearest_points_coords[residue]=nearest_points[residue].x,nearest_points[residue].y
        

    def calc_2d_forces(self,x1,y1,x2,y2,width):
        """Calculate overlap in 2D space"""
        #calculate a
        if x1>x2:
            a = x1-x2
        else:
            a = x2-x1

        a_sq=a*a
        #calculate b
        if y1>y2:
            b = y1-y2
        else: 
            b = y2-y1

        b_sq=b*b
        
        #calculate c
        from math import sqrt
        c_sq = a_sq+b_sq

        c = sqrt(c_sq)

        if c > width:
            return 0,0
        else:
            overlap = width-c
        return -overlap/2, overlap/2

    def do_step(self, values, xy_values, width=50):
        """Calculates forces between two diagrams and pushes them apart by tenth of width"""

        forces = {k:[] for k,i in enumerate(xy_values)}
        for (index1, value1), (index2,value2) in combinations(enumerate(xy_values),2):
            f = lig.calc_2d_forces(value1[0],value1[1],value2[0],value2[1],width)
            forces[index1].append(f[0])
            forces[index2].append(f[1])

        forces = {k:sum(v) for k,v in forces.items()}
        
        energy = sum([abs(x) for x in forces.values()])

        return [(forces[k]/10+v) for k, v in enumerate(values)], energy

    def make_new_projection_values(self):
        """Run do_step function until the diagramms have diverged from each other"""
        startvalues = [v for v in self.nearest_points_projection.values()]
        startvalues1 = [v for v in self.nearest_points_coords.values()]
        values = [x for x in startvalues]
        xy_values = [x for x in startvalues1]
        energy = 100
        while energy > 0.2:
            values, energy = lig.do_step(values,xy_values, width=65)
            nearest_points ={}
            nearest_points_coords={}
            xy_values=[]
            i=0
            for residue in  self.x_coordinates_from_diagramm:
                b = self.a.boundary.parallel_offset(self.closest_atoms[residue][1]*32.0+24,"left",join_style=2).convex_hull
                nearest_points[residue] = b.boundary.interpolate(values[i] % b.boundary.length)
                results_tuple = nearest_points[residue].x, nearest_points[residue].y
                xy_values.append(results_tuple)
                i+=1
        i=0
        self.new_coords_dict={}
        for residue in self.x_coordinates_from_diagramm:
            self.new_coords_dict[residue] = xy_values[i]
            i+=1




    def find_solvent_exp_atoms(self, grofile, ligand_name):
        solvent_exposed_atoms = []
        sol_exp = MDAnalysis.Universe(grofile)
        solvent_exposed = sol_exp.select_atoms("resname "+ligand_name+" and not name H* and not around 3.5 protein ")
        for atom in solvent_exposed.atoms.names:
            solvent_exposed_atoms.append(atom)
        ligand = MDAnalysis.Universe(ligand_name+".pdb")
        atom_indices=[]
        for atom in ligand.atoms.indices:
            if ligand.atoms.names[atom] in solvent_exposed_atoms:
                atom_indices.append(ligand.atoms.indices[atom])
        return atom_indices
    
    def delete_svg(self):
        list_of_plotted_svgs = []
        for residue in self.list_of_plotted_res:
            list_of_plotted_svgs.append(str(residue)+".svg")
        return list_of_plotted_svgs

    def draw_figure(self,diagramm_type, ligand_name, output_name, grofile, xtcfile):
        """Draws a clean version of the ligand (might need some """
        #Nicer ligand drawing (still needs work)
        AllChem.Compute2DCoords(self.ligand_rdkit)
        drawer = rdMolDraw2D.MolDraw2DSVG(600,300)
        opts = drawer.drawOptions()
        #Check whether solvent exposure drawing needed - at the moment just for single grofiles
        if xtcfile=="None":
            atom_indices = lig.find_solvent_exp_atoms(grofile, ligand_name)
            drawer.DrawMolecule(self.ligand_rdkit, highlightAtoms=atom_indices)
        else:
            drawer.DrawMolecule(self.ligand_rdkit)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:','')

        #write the svg figure to a file
        filesvg = open("molecule_final.svg", "w+")
        filesvg.write(svg)
        filesvg.close()

        #Change the solvent exposure colour from red to blue
        for i, line in enumerate(fileinput.input('molecule_final.svg', inplace=1)):
            sys.stdout.write(line.replace("FF7F7F","deebf7"))

        #Start by adding the bigger box
        filesvg = open("molecule_final.svg", "r")
        filesvg2 = open("molecule_modified1.svg", "w")
        lines = filesvg.readlines()
        bigger_box ="width='800px' height='600px' > <rect style='opacity:1.0;fill:white;stroke:none' width='800px' height='600px' x='0' y='0'> </rect><g transform='translate(100,150)'>'<rect style='opacity:1.0;fill:#ffffff;stroke:none' width='600' height='300' x='0' y='0'> </rect>\n"
        total= lines[:5]+[bigger_box]+lines[8:]
        filesvg2.writelines(total)
        filesvg2.close()
        filesvg.close()

        diagram="</g>"
        for residue in self.list_of_plotted_res:
            for i, line in enumerate(fileinput.input(str(residue)+'.svg', inplace=1)):
                if i <= 8:
                    continue
                else:
                    sys.stdout.write(line.replace ("</svg>","</g>"))
            for i, line in enumerate(fileinput.input(str(residue)+'.svg', inplace=1)):        
                sys.stdout.write(line.replace("</defs>","<g transform='translate("+str(int(self.new_coords_dict[residue][0])-40)+","+str(int(self.new_coords_dict[residue][1])-40)+")'>"))
            for i, line in enumerate(fileinput.input(str(residue)+'.svg', inplace=1)):        
                sys.stdout.write(line.replace("font-style:normal;","font-style:normal;font-weight:bold;"))
            for i, line in enumerate(fileinput.input(str(residue)+'.svg', inplace=1)):
                diagram = diagram+line

        #Drawing lines between res and ligand atoms (somewhat buggy)

        draw_lines = ""
        for residue in self.x_coordinates_from_diagramm:
            draw_lines=draw_lines+"<line x1='"+str(int(self.new_coords_dict[residue][0]))+"' y1='"+str(int(self.new_coords_dict[residue][1])+0)+"' x2='"+str(float(self.x_coordinates_from_diagramm[residue]))+"' y2='"+str(float(self.y_coordinates_from_diagramm[residue]))+"' style='stroke:rgb(255,0,0);stroke-width:2' />"

        filesvg = open("molecule_modified1.svg", "r")
        filesvg2 = open(output_name+".svg", "w")
        lines = filesvg.readlines()
        total1= lines[:len(lines)-1]+[diagram]+lines[len(lines)-1:]

        filesvg2.writelines(total1)
        filesvg2.close()
        filesvg.close()
        print "Drawing the figure..."


lig = LigandInteractions()

