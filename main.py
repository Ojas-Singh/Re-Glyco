# AlphaFold Test on All human proteins where Uniprot has glycosylation information.
from tqdm import tqdm
import requests
import sys, time
import os,json
import time
import config
from lib import pdb,algo,prep


default_glycan = {
    "C-linked (Man) tryptophan": "Man",
    "O-linked (Fuc...) serine": "Fuc",
    "O-linked (Fuc...) threonine": "Fuc",
    "O-linked (GalNAc...) serine": "Neu5Ac(a2-3)Gal(b1-3)GalNAc",
    "O-linked (GalNAc...) threonine": "Neu5Ac(a2-3)Gal(b1-3)GalNAc",
    "O-linked (Glc...) serine": "Glc",
    "O-linked (Glc...) threonine": "Glc",
    "O-linked (Man...) serine": "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man",
    "O-linked (Man...) threonine": "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man",
    "O-linked (Xyl...) serine": "Xyl",
    "O-linked (Xyl...) threonine": "Xyl",
    "O-linked (Xyl...) (chondroitin sulfate) serine": "GalNAc(b1-4)GlcA(b1-3)GalNAc(b1-4)GlcA(b1-3)GalNAc(b1-4)GlcA(b1-3)GalNAc(b1-4)GlcA(b1-3)Gal(b1-3)Gal(b1-4)Xyl",
    "O-linked (Xyl...) (chondroitin sulfate) threonine": "GalNAc(b1-4)GlcA(b1-3)GalNAc(b1-4)GlcA(b1-3)GalNAc(b1-4)GlcA(b1-3)GalNAc(b1-4)GlcA(b1-3)Gal(b1-3)Gal(b1-4)Xyl",
    "O-linked (GlcNAc) serine": "GlcNAc",
    "O-linked (GlcNAc) threonine": "GlcNAc",
    "N-linked (GlcNAc...) asparagine": "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    "N-linked (GlcNAc...) (complex) asparagine": "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
    "N-linked (GlcNAc...) (hybrid) asparagine": "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    "N-linked (GlcNAc...) (high mannose) asparagine": "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
}


def one_uniprot(uniprotID):
    downloadLocation = config.upload_dir
    outputFileName = uniprotID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
    requestURLpdb = f"https://alphafold.ebi.ac.uk/files/AF-{uniprotID}-F1-model_v4.pdb"
    query = requests.get(requestURLpdb, allow_redirects=True)
    outputLines = []
    downloadedLines = query.iter_lines()
    for line in downloadedLines:
        decodedLine = line.decode("utf-8")
        if decodedLine[:5] != "MODEL":
            outputLines.append(decodedLine)
    with open(outputFilePath, "w") as file:
        file.writelines("%s\n" % l for l in outputLines)
    protein = pdb.parse(outputFilePath)
    df = pdb.to_DF(protein)
    uniprotRequestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprotID}"
    uniprotResponse = requests.get(
        uniprotRequestURL, headers={"Accept": "application/json"}
    )
    if not uniprotResponse.ok:
        uniprotResponse.raise_for_status()
        sys.exit()
    uniprotResponseJSON = uniprotResponse.json()
    uniprotFeatures = uniprotResponseJSON["features"]
    selectedGlycans = {}
    uniprot_CARBOHYD = {}
    for item in uniprotFeatures:
        if item["type"] == "CARBOHYD":
            uniprot_CARBOHYD=item
            resname = df.loc[df['ResId'] == int(item["begin"]), 'ResName'].iloc[0]
            key = str(int(item["begin"])) + "_A"
            info = item["description"]
            
            # Check if the info matches any key in the glycan_data and get the corresponding value
            value = default_glycan.get(info, None)
            
            if value:  # If value is not None, add to the selectedGlycans dictionary
                selectedGlycans[key] = value
    # Extract glycans and glycosylation_locations from selectedGlycans
    glycans = [glycan for glycan in selectedGlycans.values() if glycan]  # Filter out empty glycans
    glycosylation_locations = [location for location in selectedGlycans.keys()]

    # Assuming you have the necessary pdb and algo modules and methods
    protein = pdb.parse(f'{downloadLocation}{uniprotID}.pdb')
    g, clash, box = algo.attach(protein, glycans, glycosylation_locations)
    outputshortfilepath = f'{uniprotID}_glycosylated.pdb'
    outputfilepath2 = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.exportPDB(outputfilepath2, pdb.to_normal(g))
    return clash, box ,uniprot_CARBOHYD


def process_all_uniprots(file_path):
    results = []
    
    # Load all Uniprot IDs into a list
    with open(file_path, 'r') as file:
        uniprot_ids = [line.strip() for line in file]
    
    # Process each Uniprot ID with a progress bar
    for uniprotID in tqdm(uniprot_ids, desc="Processing Uniprot IDs"):
        try:
            clash_data, box_data, uniprot_CARBOHYD = one_uniprot(uniprotID)
            results.append({
                'uniprotID': uniprotID,
                'clash_data': clash_data,
                'box_data': box_data,
                'uniprot_CARBOHYD': uniprot_CARBOHYD
            })
        except Exception as e:
            print(f"Error processing {uniprotID}: {e}")

    # Ensure the output directory exists, if not, create it
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    # Save the results to an output file in the specified directory
    output_file_path = os.path.join(OUTPUT_DIR, 'processed_results.json')
    with open(output_file_path, 'w') as outfile:
        json.dump(results, outfile)

    return results



if __name__ == "__main__":
    OUTPUT_DIR = 'results/'
    file_path = 'human_glycoproteins_test.txt'
    results = process_all_uniprots(file_path)
    print("Processing complete.")