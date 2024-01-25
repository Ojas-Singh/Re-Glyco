import numpy as np
import MDAnalysis as mda
import shlex
import os,subprocess,glob
from shutil import copyfile,move

AMINO_SUBSTITUTION = {'AIB': 'ALA', 'ALM': 'ALA', 'AYA': 'ALA', 'BNN': 'ALA', 'CHG': 'ALA', 'CSD': 'ALA', 'DAL': 'ALA', 'DHA': 'ALA', 'DNP': 'ALA', 'FLA': 'ALA', 'HAC': 'ALA', 'MAA': 'ALA', 'PRR': 'ALA', 'TIH': 'ALA', 'TPQ': 'ALA', 'AGM': 'ARG', 'DAR': 'ARG', 'HAR': 'ARG', 'MMO': 'ARG', 'ARM': 'ARG', 'ARN': 'ARG', 'HMR': 'ARG', 'ACL': 'ARG', 'MEN': 'ASN', 'DSG': 'ASN', 'DSP': 'ASP', 'BHD': 'ASP', '2AS': 'ASP', 'ASQ': 'ASP', 'ASB': 'ASP', 'ASA': 'ASP', 'ASK': 'ASP', 'ASH': 'ASP', 'ASL': 'ASP', 'DAS': 'ASP', 'BCS': 'CYS', 'BUC': 'CYS', 'C5C': 'CYS', 'C6C': 'CYS', 'CCS': 'CYS', 'CEA': 'CYS', 'CME': 'CYS', 'CSO': 'CYS', 'CSP': 'CYS', 'CSS': 'CYS', 'CSW': 'CYS', 'CSX': 'CYS', 'CY1': 'CYS', 'CY3': 'CYS', 'CYG': 'CYS', 'CYM': 'CYS', 'CYP': 'CYS', 'CYQ': 'CYS', 'CYX': 'CYS', 'DCY': 'CYS', 'EFC': 'CYS', 'OCS': 'CYS', 'PEC': 'CYS', 'PR3': 'CYS', 'SCH': 'CYS', 'SCS': 'CYS', 'SCY': 'CYS', 'SHC': 'CYS', 'SMC': 'CYS', 'SOC': 'CYS', 'GLH': 'GLU', 'GGL': 'GLU', 'PCA': 'GLU', '5HP': 'GLU', 'DGL': 'GLU', 'CGU': 'GLU', 'GMA': 'GLU', ('DGN', 'X'): 'GLN', 'GLZ': 'GLY', 'SAR': 'GLY', 'NMC': 'GLY', 'GL3': 'GLY', 'GSC': 'GLY', 'MPQ': 'GLY', 'MSA': 'GLY', 'DHI': 'HIS', 'HID': 'HIS', 'HIC': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS', 'MHS': 'HIS', 'NEM': 'HIS', 'NEP': 'HIS', '3AH': 'HIS', 'DIL': 'ILE', 'IIL': 'ILE', 'BUG': 'LEU', 'NLE': 'LEU', 'NLP': 'LEU', 'NLN': 'LEU', 'DLE': 'LEU', 'CLE': 'LEU', 'MLE': 'LEU', 'LYM': 'LYS', 'ALY': 'LYS', 'LYZ': 'LYS', 'LYN': 'LYS', 'LLY': 'LYS', 'LLP': 'LYS', 'SHR': 'LYS', 'TRG': 'LYS', 'DLY': 'LYS', 'KCX': 'LYS', 'FME': 'MET', 'CXM': 'MET', 'OMT': 'MET', 'MSE': 'MET', 'DAH': 'PHE', 'DPN': 'PHE', 'HPQ': 'PHE', 'PHI': 'PHE', 'PHL': 'PHE', 'DPR': 'PRO', 'HYP': 'PRO', 'OAS': 'SER', 'MIS': 'SER', 'SAC': 'SER', 'SVA': 'SER', 'SET': 'SER', 'SEP': 'SER', 'SEL': 'SER', 'DSN': 'SER', 'ALO': 'THR', 'BMT': 'THR', 'DTH': 'THR', 'THO': 'THR', 'TPO': 'THR', 'DTR': 'TRP', 'HTR': 'TRP', 'LTR': 'TRP', 'TPL': 'TRP', 'TRO': 'TRP', 'DTY': 'TYR', 'IYR': 'TYR', 'PAQ': 'TYR', 'PTR': 'TYR', 'STY': 'TYR', 'TYB': 'TYR', 'TYM': 'TYR', 'TYO': 'TYR', 'TYQ': 'TYR', 'TYS': 'TYR', 'TYY': 'TYR', 'DIV': 'VAL', 'DVA': 'VAL', 'MVA': 'VAL'}

def calc_SASA(pdbname,xtcname,index,groupname,selectgroups,outfile,outfileres,outfileatom,probe,ndots,endtime):

   cmd = f"gmx sasa -f {xtcname} -s {pdbname} -n {index} -o {outfile} -or {outfileres} " \
         f"-oa {outfileatom} -surface '{groupname}' -xvg none -probe {probe} " \
         f"-output '{selectgroups}' -ndots {ndots} -e {endtime}"
   process = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
   
   if process.returncode != 0:
      print(process.stderr)
      raise SystemExit("Error in Gromacs, stopping...")

   # Read in the Gromacs output files
   sasas = np.loadtxt(outfileres)
   sasaa = np.loadtxt(outfileatom)
   return  sasas, sasaa, process.returncode


   
def run(pdbfiles,xtcfiles,probe,ndots,endtime,outfile):

      outrelativesasa=[]
      outrelativesasaaa=[]      

      for iglycan in range(len(xtcfiles)):
          
          pdb=pdbfiles[iglycan]
          xtc=xtcfiles[iglycan]
          u=mda.Universe(pdb,xtc)
          
          tmpindex='index0.ndx'
          tmpsel='prot; gly'
          tmpsys='system'
          baresel='prot'
          tmpsasa  = 'sasa1.xvg'
          tmpsasar = 'sasar.xvg'
          tmpsasaa = 'sasaa.xvg'
          tmppdb = 'test1.pdb'
          tmpbarepdb = 'test2.pdb'
          tmpbarextc = 'test2.xtc'
          with open(tmppdb,'w') as f:
             with open(pdb,'r') as g:
                for line in g: 
                   if "ATOM" in line:
                      myres = line[17:20]
                      if myres in AMINO_SUBSTITUTION.keys():
                        line.replace(myres,AMINO_SUBSTITUTION[myres])
                      f.write(line)
          u=mda.Universe(tmppdb,xtc)
          sel_P=u.select_atoms('protein')
          sel_P.write(tmpbarepdb)
          with open(tmpbarepdb) as fff:
            with open('ggg','w') as ggg:
               for line in fff:
                  if "ATOM" in line:
                     ggg.write(line)    
          move('ggg',tmpbarepdb)
          with mda.coordinates.XTC.XTCWriter(tmpbarextc,n_atoms=sel_P.atoms.n_atoms) as w:
               for tp in u.trajectory:
                  w.write(sel_P.atoms)
          sel_G=u.select_atoms('not protein')
          with mda.selections.gromacs.SelectionWriter(tmpindex) as w:
              w.write(sel_P,name='prot')
              w.write(sel_G,name='gly')
              w.write(u.atoms,name='system')
          sasar,sasaa,_ = calc_SASA(tmppdb,xtc,tmpindex,tmpsys,tmpsel,tmpsasa,tmpsasar,tmpsasaa,probe,ndots,endtime)
          if iglycan==0:
              baresasar,baresasaa,_ = calc_SASA(tmpbarepdb,tmpbarextc,tmpindex,baresel,baresel,tmpsasa,tmpsasar,tmpsasaa,probe,ndots,endtime)
          baresasar_idx=np.where(baresasar[:,1]!=0)
          occupancy_r=np.zeros(len(baresasar[:,1]))
          occupancy_r[baresasar_idx]=1.0
          relativesasa=np.zeros(len(baresasar[:,1]))
          relativesasa[baresasar_idx]=(baresasar[:,1][baresasar_idx]-sasar[:baresasar.shape[0],1][baresasar_idx])/baresasar[:,1][baresasar_idx]
          baresasaa_idx=np.where(baresasaa[:,1]!=0)          
          occupancy_a=np.zeros(len(baresasaa[:,1]))
          occupancy_a[baresasaa_idx]=1.0
          relativesasaaa=np.zeros(len(baresasaa[:,1]))
          relativesasaaa[baresasaa_idx]=(baresasaa[:,1][baresasaa_idx]-sasaa[:baresasaa.shape[0],1][baresasaa_idx])/baresasaa[:,1][baresasaa_idx]
          outrelativesasa.append(relativesasa)
          outrelativesasaaa.append(relativesasaaa)
          os.remove(tmpindex)
          os.remove(tmpsasa)
          os.remove(tmpsasar)
          os.remove(tmpsasaa)
          for afile in glob.glob("#sasa*"):
             os.remove(afile)
      outrelativesasa=np.array(outrelativesasa)
      outrelativesasaaa=np.array(outrelativesasaaa)
      meanSASA=np.mean(outrelativesasa,axis=0)
      unw=mda.Universe(tmppdb)
      sel_Pnw=unw.select_atoms('protein')
      for iresidue in range(len(sel_P.residues)):
            sel_Pnw.residues[iresidue].atoms.tempfactors=meanSASA[iresidue]*100
            sel_Pnw.residues[iresidue].atoms.occupancies=occupancy_r[iresidue]
      sel_Pnw.write(outfile.format(probe))
      