#########################################################################################################################################
#########################################################################################################################################
############ This code was originally used as a module for the 4th stage of my protein-docking pipeline, stages 1-3 being
############ 1) selection of structures and prep, 2) docking of structures, 3) clustering of docking output ("decoys"), which
############ generally consisted of 10,000 structures and anywhere from 1 to 20 clusters. The 4th stage was, for each cluster,
############ to select a "representative structure" on which to perform additional analysis. An optional extra part of this representative structure
############ selection was to extract from the rep a consensus set of residues ("core"), or those that fell within a certain distance threshold
############ compared to other members of the cluster. When both the rep structure and its "core" are loaded into a molecular visualization program,
############ the regions of the protein most variable across decoys, which often correspond to regions most active during docking, can be clearly seen
############ as those portions of the rep structure NOT present in the core.
############ The file "huLEDGF_vs_hivINTEGRASEmutant6.pdb" and "huLEDGF_vs_hivINTEGRASEmutant6_core.pdb" show this. The former is included in my email,
############ while the latter is generated when the code below is executed 
#########################################################################################################################################
#########################################################################################################################################


import os
import sys
##### Define working directory
wdir='C:\\Users\\Norman\\Documents\\docking'
pdbdir=wdir+'\\'+'example_docking_structures2'
##### Bio.PDB is part of the Biopython package, so it should load if Biopython is installed 
import Bio.PDB.PDBParser
parser=Bio.PDB.PDBParser()

class docking_ensemble():
    ####################################################################################################################################
    ########## Loop imports pdb structures.
    ##### Input is a list of pdb structures to be analyzed.
    ##### Output is a dictionary (pdbstructures) containing Bio.PDB structure objects as values, with pdb filenames (from input) as keys
    def __init__(self,runid,pdbdir,pdbfilenames):
        self.runid=runid
        pdbstructures=dict()
        for i,pdbfilename in enumerate(pdbfilenames):
            print pdbfilename
            pdbname=pdbfilename.split('.')[0]
            structure=parser.get_structure(pdbname,pdbdir+'\\'+pdbfilename)
            pdbstructures[pdbfilename.split('.')[0]]=structure
        self.sraw=pdbstructures
    ####################################################################################################################################
    ########## The raw Bio.PDB data object (self.sraw) is an iterably-callable tree-structure
    ########## where you go first from structure, to chain to residue, to atom, using a .children() method function... not very convenient
    ########## if you want to go in and grab residues on an as-needed basis. This function converts that Bio.PDB data object to a more flat structure
    ##### Input is a dict of Bio.PDB structures (self.sraw) and a pdb:chain dict (param2)... Ex. param 2: {'1KFX':['L','S'], '1KFU':['L','S']}
    ##### Output is a dictionary that is "residue centric" (self.sdict), with keys that are of the form: *chain*_*atomnumber* 
    def struct_toresidues(self,pdbchains):
        pdbresdict=dict()
        self.chains=pdbchains
        for pdbname in self.sraw.keys():
            structure=self.sraw[pdbname]
            # If an X-Ray or HADDOCK docked structure, index will always be 0, but if NMR, may be greater
            # This code doesn't deal with ensembles of ensembles, so always 0
            model=structure[0]
            resdict=dict()
            chainids=self.chains[pdbname]
            ## Loops that collapses chain and residue names into one key, stores atoms as dict
            for chainid in chainids:
                chain=model[chainid]
                for res in chain:
                    resdict[chainid+'_'+str(res.get_id()[1])]=res
            pdbresdict[pdbname]=resdict
        self.sdict=pdbresdict
    ####################################################################################################################################
    ########## Calculates the structure of the ensemble with the lowest cumulative rmsd vs other structures
    ########## Stores all pairwise rmsd values
    ##### Input is self.sdict
    ##### Output is most central structure (self.best) and all pairwise rmsd values (self.rmsd)
    def find_representative_pdb(self,chainmap):
        rmsd=dict()
        # Initialize minimum rmsd to very high value
        bestrmsd=float(1000)
        for pdbname1 in self.sdict.keys():
            # rmsdsum1 stores cumulative rms for pdbname1
            rmsdsum1=float(0)
            rmsd1=dict()
            for pdbname2 in self.sdict.keys():
                if pdbname2==pdbname1:
                    continue
                chainmap2=chainmap[pdbname1][pdbname2]
                fixed,moving=get_fixed_moving_CAlists(self.sdict[pdbname1],self.sdict[pdbname2],chainmap2,0,False)
                try:
                    rms=get_rms(fixed,moving)
                except ZeroDivisionError:
                    rms=float(0)
                rmsd1[pdbname2]=rms
                rmsdsum1+=rms
            rmsd[pdbname1]=rmsd1
            print pdbname1+' has cumulative RMSD: '+str(rmsdsum1)
            # Promotes pdbname1 to bestpdb if cumulative rmsd from rmsd_sum is new minimum
            if rmsdsum1 < bestrmsd:
                bestpdb=pdbname1
                bestrmsd=rmsdsum1
        self.chainmap=chainmap
        self.best=bestpdb
        self.rmsd=rmsd

class align_ensemble(docking_ensemble):
    ####################################################################################################################################
    ########## Superimposes all structures to the "best" or most central structure, self.best
    ##### Input is self.sraw, self.sdict, self.best, and an output filename
    ##### Output is a set of new pdb files, one for each structure in the docking ensemble aligned to the most central
    ##### Output does not save aligned structures as a variable, since the output files could be used to initialize a new instance of the class
    def __init__(self,runid,pdbdir,pdbfilenames,pdbchains,chainmap):
        docking_ensemble.__init__(self,runid,pdbdir,pdbfilenames)
        self.struct_toresidues(pdbchains)
        self.find_representative_pdb(chainmap)
    def align_all_to_common_structure(self):
        io=Bio.PDB.PDBIO()
        superimposer=Bio.PDB.Superimposer()
        for pdbname in self.sdict.keys():
            if pdbname==self.best:
                continue
            print pdbname
            chainmap2=self.chainmap[self.best][pdbname]
            fixed,moving=get_fixed_moving_CAlists(self.sdict[self.best],self.sdict[pdbname],chainmap2,0,False)
            try:
                superimposer.set_atoms(fixed,moving)
                moving_structure=self.sraw[pdbname]
                superimposer.apply(moving_structure.get_atoms())
            except ZeroDivisionError:
                a=1
            outfile=open(pdbdir+'\\'+pdbname+'_aligned.pdb','w')
            io.set_structure(moving_structure)
            io.save(outfile)
            outfile.close()
    ####################################################################################################################################
    ########## Extracts, from the best structure of the ensemble, the "best" residues - those that are within a certain average distance threshold from the remainder of the ensemble
    ##### Input is bestpdb (not self.bestpdb because bestpdb comes from the original, unaligned docking_ensembl object), a distance cutoff, self.sraw, self.sraw, self.chains,
    ##### Output is a new .pdb file with only the core of "best" residues from the bestpdb
    def extract_core(self,pdbdir,bestpdb,distcutoff):
        io=Bio.PDB.PDBIO()
        beststruct=self.sraw[bestpdb]
        model=beststruct[0]
        bestres=self.sdict[bestpdb]
        # Iterate through residues in beststruct and calculate by-residue cumulative rmsd
        distres=dict()
        for pdbname in self.sdict.keys():
            if pdbname==bestpdb:
                continue
            chainmap2=self.chainmap[self.best][pdbname]
            fixed,moving=get_fixed_moving_CAlists(self.sdict[self.best],self.sdict[pdbname],chainmap2,0,True)
            for key in fixed.keys():
                atom1_coord=fixed[key].get_coord()
                atom2_coord=moving[key].get_coord()
                dist=euclidean_distance(atom1_coord,atom2_coord)
                try:
                    dists=distres[key]
                except KeyError:
                    dists=[]
                dists.append(dist)
                distres[key]=dists
        # Iterate through residues in beststruct and remove those not, on average, within the cutoff distance
        for key in distres.keys():
            dists=distres[key]
            avgdist=average(dists)
            if avgdist > distcutoff:
                chainid=key.split('_')[0]
                pos=key.split('_')[1]
                chain=model[chainid]
                res=bestres[key]
                chain.detach_child(res.id)
        outfile=open(pdbdir+'\\'+bestpdb+'_core.pdb','w')
        io.set_structure(beststruct)
        io.save(pdbdir+'\\'+bestpdb+'_core.pdb')
        outfile.close()

########## Grabs pairs of alpha-carbon chains to be aligned
##### Input is pairs of resdicts (as generated by make_residues_dict) as params 1 and 2, a chainmap saying which chain in dict1 correspond to which chain in
##### ... chain 2 as param 3, and on offset value in case the residue numbering system in chain2 is ahead (+) or behind (-) compared to chain1
##### Output is two lists of alpha carbon atoms, the first being the reference and the second being that which is to be superimposed over reference
def get_fixed_moving_CAlists(resdict1,resdict2,chainmap,offset,asdict):
    fixed=[]
    moving=[]
    if asdict:
        fixed=dict()
        moving=dict()
    for key in resdict1.keys():
        res1=resdict1[key]
        try:
            chain1,pos1=key.strip().split('_')
            chain2=chainmap[chain1]
            pos2=str(int(pos1)+offset)
            key2=chain2+'_'+pos2
            res2=resdict2[key2]
        # KeyError will typically occur when one structure contains a residue that the other doesn't as a result of being slightly longer
        except KeyError:
            continue
        if not asdict:
            try:
                fixed.append(res1['CA'])
                moving.append(res2['CA'])
            # KeyError will typically occur of the residue has incomplete X-ray or NMR data, and thus lacks an alpha carbon
            except KeyError:
                continue
        else:
            try:
                fixed[key]=res1['CA']
                moving[key]=res2['CA']
            except KeyError:
                continue
    return [fixed,moving]
    
##### Unfortunately rmsd calculation can only be done pairs at a time (as far as I know),
##### so here, submit a pair of residue dictionaries, each corresponding to a
##### pdb file, at a time.
def get_rms(fixed,moving):
    superimposer=Bio.PDB.Superimposer()
    superimposer.set_atoms(fixed,moving)
    return superimposer.rms

def euclidean_distance(atom1_coord,atom2_coord):
    x1,y1,z1=atom1_coord
    x2,y2,z2=atom2_coord
    distance=((x2-x1)**2.0+(y2-y1)**2.0+(z2-z1)**2.0)**0.5
    return distance

def average(numbers):
    return sum(numbers)/float(len(numbers))

##### Initializes data used as input to create a docking_ensemble or aligned_ensemble object instance
def initialize_dockingrun_information(runid,decoys,aligned):
    if not aligned:
        docking_runid_tags=[runid+str(n) for n in decoys]
    else:
        docking_runid_tags=[runid+str(n)+'_aligned' for n in decoys]
    docking_filenames=[drtag+'.pdb' for drtag in docking_runid_tags]
    ##### Grow a chainmap dictionary where all chains 'A' correspond to each other across all structures, and all chains 'B' likewise
    ##### In a different scenario, i.e. an ensembl of docked structures from different source .pdb files, chains may be different
    chains=dict(zip(docking_runid_tags,['A','B']*len(decoys)))
    chainmap=dict()
    for dridtag1 in docking_runid_tags:
        chainmap1=dict()
        for dridtag2 in docking_runid_tags:
            if dridtag2==dridtag1:
                continue
            chainmap2=dict()
            chainmap2['A']='A'
            chainmap2['B']='B'
            chainmap1[dridtag2]=chainmap2
        chainmap[dridtag1]=chainmap1
    return [docking_filenames,chains,chainmap]

#####################################################################################################################################
#####  PRE-ALIGNMENT (determine representative structure)
runid='huLEDGF_vs_hivINTEGRASEmutant'
decoys=range(1,10)
aligned=False
docking_filenames,chains,chainmap=initialize_dockingrun_information(runid,decoys,aligned)
print "\nCreate align_ensemble instance with docking information (ensembl.best contains most representative docking decoy) for "+str(len(decoys))+" unaligned decoys"
ensembl_aln=align_ensemble(runid,pdbdir,docking_filenames,chains,chainmap)
print "\nAligning "+str(len(docking_filenames))+" docking decoys"
ensembl_aln.align_all_to_common_structure()

#####################################################################################################################################
#####  POST-ALIGNMENT (align all to representative, extract "core" of consensus residues from representative
aligned=False
docking_filenames,chains,chainmap=initialize_dockingrun_information(runid,decoys,aligned)
coredistcutoff=float(0.02)
print "\nCreate align_ensemble instance with docking information for "+str(len(decoys))+" decoys aligned to "+ensembl_aln.best
ensembl_core=align_ensemble(runid,pdbdir,docking_filenames,chains,chainmap)
print "\nExtracting core of < "+str(coredistcutoff)+" angstroms consensus from representative structure: "+ensembl_aln.best
ensembl_core.extract_core(pdbdir,ensembl_aln.best,coredistcutoff)








################# creates docking_ensemble object (parent of align_ensemble) - not needed ##################

####### Read in raw docking .pdb structural info   
##ensembl=docking_ensemble(runid,pdbdir,docking_filenames)
####### "Unpack" structural info from iterably-callable tree format; store residue 3D info in accessible format
##ensembl.struct_toresidues(chains)
####### Find the structure that is the best representative of the ensemble
##ensembl.find_representative_pdb(chainmap)
