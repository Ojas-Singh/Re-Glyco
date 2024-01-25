# Re-Glyco : A GlycoProtein Builder

Re-Glyco is a tool we designed to restore the missing glycosylation on glycoproteins deposited in the RCSB PDB or in the EBI-EMBL AlphaFold protein structure database. To get started, upload your protein structure file or choose a pre-existing AlphaFold or PDB structure, and let Re-Glyco do the rest!
Currently supported function includes :
- N-GlcNAcylation
- O-GalNAcylation
- O-GlcNAcylation
- O-Fucosylation
- O-Mannosylation
- O-Glucosylation
- O-Xylosylation
- C-Mannosylation

This tool is currently hosted under GlycoShape project and can be accessed at https://glycoshape.org/reglyco

# Installation

```
sudo apt install build-essential
conda create -n reglyco python=3.12
conda activate reglyco
conda install conda-forge::gromacs 
pip install -r requirements.txt

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup update

cd glycors
maturin develop --release

```
# Running the API 

Modify the config.py to locate the GlycoShape Database directory.

```
gunicorn -w 4 api:app --timeout 900
gunicorn -w 4 api:app --timeout 1500 -b 127.0.0.1:8001
```

# Citation

All of the data provided is freely available for academic use under Creative Commons Attribution 4.0 (CC BY-NC-ND 4.0 Deed) licence terms. Please contact us at elisa.fadda@mu.ie for Commercial licence. If you use this resource, please cite the following papers:

Callum M Ives and Ojas Singh et al. Restoring Protein Glycosylation with GlycoShape bioRxiv (2023).  https://doi.org/10.1101/2023.12.11.571101

# Future roadmap
- CLI interface
- Making everything fast.




