# Re-Glyco

The Python software utilizes molecular dynamics (MD) simulation results from the glycoshape database to re-glycosylate protein structures predicted by alphafold, a deep learning method for protein structure prediction. The input to the software is a protein structure output by alphafold, and the output is a modified protein structure with glycans (sugar molecules) added at appropriate sites. The incorporation of glycans is achieved through the use of MD simulation results, which ensures that the resulting glycosylated protein structure is physically realistic and stable. This software can be useful for researchers studying the role of glycans in protein function and the impact of glycosylation on protein structure and stability.


# Installation
```
conda install -n reglyco python=3.10
conda activate reglyco
python -m pip install -r requirements.txt
streamlit run main.py
```


To implement!
1. Custom Glycosylation spot.
2. Clash Detector warning.
2. Advanced method.
    wiggle[need work] and different cluster search.

3. best way is to find psi phi first then use that value with wiggle if thats fails too wiggle low confidence area.

Automate Torsion finder.
