# Re-Glyco : A GlycoProtein Builder

[Re-Glyco](https://glycoshape.org/reglyco) is a tool we designed to restore the missing glycosylation on glycoproteins deposited in the RCSB PDB or in the EBI-EMBL AlphaFold protein structure database. To get started, upload your protein structure file or choose a pre-existing AlphaFold or PDB structure, and let Re-Glyco do the rest!
Currently supported function includes :
- N-GlcNAcylation
- O-GalNAcylation
- O-GlcNAcylation
- O-Fucosylation
- O-Mannosylation
- O-Glucosylation
- O-Xylosylation
- C-Mannosylation

This tool is currently hosted under [GlycoShape project](https://glycoshape.org/)

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
```

if you want to run on different port
```
gunicorn -w 4 api:app --timeout 1500 -b 127.0.0.1:8001 
```

# Frontend

The current version of Re-Glyco only has backend API service which can be used by [GlycoShape Website](https://glycoshape.org) code of the frontend website is [here](https://github.com/Ojas-Singh/GlycoShape).

You can also use python script to access the API upload your pdb file and glycosylate it, a dummy code is provided as demo.py in this repo.


# Demo
Once the API is running locally at your desired port.
change the API_BASE_URL = "https://glycoshape.org"  in demo.py to "http://127.0.0.1:8000"

and run using 
```
python demo.py
```

This will upload the [pdb](AF-P29016-F1-model_v4.pdb) to the API and will produced glycosylated structures can be found at output folder or temp_files.

It can take upto 5 minute for a job.

# Citation

All of the data provided is freely available for academic use under Creative Commons Attribution 4.0 (CC BY-NC-ND 4.0 Deed) licence terms. Please contact us at elisa.fadda@mu.ie for Commercial licence. If you use this resource, please cite the following papers:

Callum M Ives and Ojas Singh et al. Restoring Protein Glycosylation with GlycoShape [bioRxiv (2023)](https://www.biorxiv.org/content/10.1101/2023.12.11.571101v1.full).

# Future roadmap
- CLI interface
- Density fitting
- Fitness function with Non-bonded interaction.




