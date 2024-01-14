from flask import Flask, request, jsonify, make_response
from flask_cors import CORS
import requests
import sys, time
import os,json
import time
from datetime import datetime
import config
from lib import pdb,algo,prep
from thefuzz import fuzz
import re




app = Flask(__name__)

app.config['MAX_CONTENT_LENGTH'] = 16 * 1000 * 1000
CORS(app)
CORS(app, resources={r"/api/*": {"origins": "*"}})
CORS(app, supports_credentials=True)


def load_glycan_data(directory_path):
    glycan_data = []
    wurcs_with_filename = []
  
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".json"):
                with open(os.path.join(root, file), "r") as f:
                    try:
                        data = json.load(f)
                        wurcs = data.get("wurcs")  # Using get() avoids a KeyError if 'wurcs' is not present
                        if wurcs is not None:  # Check to ensure 'wurcs' was present in the data
                            wurcs_with_filename.append((wurcs, data.get("iupac")))
                            glycan_data.append((data, data.get("iupac")))  # Storing data with filename
                    except:
                        pass
    with open(directory_path+'GLYCAN_TYPE.json', 'r') as file:
        data = json.load(file)
    N_Glycans =  data.get('N', [])
    O_Glycans =  data.get('O', [])
    Oligomannose = data.get('Oligomannose',[])
    GAGs= data.get('GAG', [])
    Complex = data.get('Complex',[])
    Hybrid = data.get('Hybrid',[])
   
    return glycan_data, wurcs_with_filename, N_Glycans, O_Glycans, GAGs , Oligomannose, Complex, Hybrid


glycans_with_filename, wurcs_with_filename, N_Glycans, O_Glycans, GAGs, Oligomannose, Complex, Hybrid = load_glycan_data(config.data_dir)

def get_filtered_glycanlist(path, residue_name):
    try:
        all_folders = os.listdir(path)
    except FileNotFoundError:
        print(f"The directory {path} was not found.")
        return []
    match residue_name:
        case "ASN":
            filtered_folders = N_Glycans
            return filtered_folders
        case "SER":
            filtered_folders = [folder for folder in all_folders if folder.endswith("GalNAc") or folder.endswith("Fuc") or folder.endswith("Glc") or folder.endswith("Xyl")]
            filtered_folders.append("GlcNAc")
            return filtered_folders
        case "THR":
            filtered_folders = [folder for folder in all_folders if folder.endswith("GalNAc") or folder.endswith("Fuc") or folder.endswith("Man")]
            filtered_folders.append("GlcNAc")
            return filtered_folders
        # case "TYR":
        #     filtered_folders = ["Man"]
        #     return filtered_folders
        case "TRP":
            filtered_folders = ["Man"]
            return filtered_folders
        case _:
            return all_folders
        





def score_glycan(glycan_with_filename, search_string):
    glycan_data, filename = glycan_with_filename
    # Considering using the token_set_ratio which handles partial string matches well
    return fuzz.token_set_ratio(str(glycan_data).lower(), search_string.lower())



def search_glycans(search_string, glycans_with_filename):
    scored_glycans = [(score_glycan(glycan_with_filename, search_string), glycan_with_filename) for glycan_with_filename in glycans_with_filename]
    scored_glycans.sort(reverse=True, key=lambda x: x[0])
    return [filename for score, (glycan, filename) in scored_glycans if score > 80]  # Consider only items with a score above a threshold



@app.route('/api/available_glycans', methods=['GET'])
def list_glycans():
    base_dir = config.data_dir
    if os.path.exists(base_dir):
        dir_list = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
        return jsonify({'glycan_list': dir_list})
    else:
        return jsonify({'error': 'Directory not found'})
    
with open(f'{config.data_dir}GLYCOSHAPE.json', 'r') as file:
    GDB_data = json.load(file)

@app.route('/api/fetch_glytoucan', methods=['GET'])
def fetch_glytoucan():
    # Retrieve the ID from the request arguments
    entry_id = request.args.get('id')

    # Search for the entry with the given ID
    entry = GDB_data.get(entry_id)

    # Check if the entry exists
    if entry:
        return jsonify(entry)
    else:
        return jsonify({"error": "Entry not found"}), 404
    

@app.route('/api/search', methods=['POST'])
def search():
    string = request.json.get('search_string')
    
    # Ensure the search string is not None
    if string is None:
        return jsonify({'error': 'Search string is None'})

    # Clean up the search string
    string = string.strip()

    # Base directory path
    base_dir = '/mnt/database/DB'

    # Check if base directory exists
    if not os.path.exists(base_dir):
        return jsonify({'error': 'Directory not found'})

    # If the search string is empty, return the list of directories
    if string == "all":
        dir_list = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
        return jsonify({'search_string': string, 'results': dir_list})
    else:
        # Based on the specific search strings, return appropriate lists
        if string == "N-Glycans":
            # Your logic to get N-Glycans list
            l = N_Glycans
        elif string == "O-Glycans":
            # Your logic to get O-Glycans list
            l = O_Glycans
        elif string == "GAGs":
            # Your logic to get GAGs list
            l = GAGs
        elif string == 'Oligomannose':
            l = Oligomannose
        elif string == 'Complex':
            l = Complex
        elif string == 'Hybrid':
            l = Hybrid
        else:
            # For other search strings
            if os.path.exists(config.data_dir):
                l = search_glycans(string, glycans_with_filename)
            else:
                return jsonify({'error': 'Directory not found'})

        return jsonify({'search_string': string, 'results': l})

        

def wurcs_search(search_wurcs, wurcs_with_filename):
    scores = []
    for wurcs, filename in wurcs_with_filename:
        score = fuzz.partial_ratio(str(search_wurcs).lower(), wurcs.lower())
        scores.append((score, filename))
    scores.sort(reverse=True, key=lambda x: x[0])
    top_filenames = [filename for score, filename in scores[:10]]
    return top_filenames

@app.route('/api/wurcs', methods=['POST'])
def wurcs():
    string = request.json.get('wurcs_string')

    if os.path.exists( config.data_dir):
        l = wurcs_search(string, wurcs_with_filename)
        return jsonify({'search_string':string ,'results': l})
    else:
        return jsonify({'error': 'Directory not found'})


@app.route('/api/rcsb', methods=['POST'])
def rcsb():
    downloadLocation = config.upload_dir
    pdbID = request.json.get('uniprot')
    requestURLcif = f"https://files.rcsb.org/download/{pdbID}.cif"
    outputFileName = pdbID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
    requestURLpdb = f"https://files.rcsb.org/download/{pdbID}.pdb"
    response = requests.get(requestURLpdb)
    
    if response.status_code == 200:
        with open(f"{downloadLocation}{pdbID}.pdb", "wb") as f:
            f.write(response.content)
        print(f"Downloaded {pdbID}.pdb successfully!")
        protein = pdb.parse(outputFilePath)
        df = pdb.to_DF(protein)
        final = []
        # sequence, shift = pdb.get_sequence_from_pdb(filename)
        # spots = pdb.find_glycosylation_spots(sequence,shift)

        sequence_with_info, sequences = pdb.get_sequence_from_pdb(outputFilePath )
        spots = pdb.find_glycosylation_spots(sequence_with_info)
        for i in spots:
                resname = df.loc[(df['ResId'] == int(i[1])) & (df['Chain'] == i[2]), 'ResName'].iloc[0]
                output = {  
                        'residueTag': int(i[0]),
                        'residueID': int(i[1]),
                        'residueName': resname,
                        'residueChain': i[2],
                        'glycanIDs':  get_filtered_glycanlist(config.data_dir,resname),
                            }  
                final.append(output)
        data = {
            "sequenceLength": len(sequence_with_info),
            "sequence": sequences ,
            "glycosylations": [],
        }
        return jsonify({'uniprot':f'{pdbID}','glycosylation_locations': data, 'requestURL': requestURLcif, 'configuration' : final})

    else:
        print(f"Failed to download {pdbID}.pdb. HTTP Status Code: {response.status_code}")
        response.raise_for_status()
        sys.exit()
    
    

    
@app.route('/api/uniprot', methods=['POST'])
def uniprot():
    # base_dir = config.data_dir
    # if os.path.exists(base_dir):
    #     dir_list = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    downloadLocation = config.upload_dir
    uniprotID = request.json.get('uniprot')
    requestURLcif = f"https://alphafold.ebi.ac.uk/files/AF-{uniprotID}-F1-model_v4.cif"
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
    uniprotGlycosylations = []
    uniprotSequence = uniprotResponseJSON["sequence"]
    uniprotFeatures = uniprotResponseJSON["features"]
    uniprot_spots= []
    final = []
    for item in uniprotFeatures:
        if item["type"] == "CARBOHYD":
            uniprotGlycosylations.append(item)
            resname = df.loc[df['ResId'] == int(item["begin"]), 'ResName'].iloc[0]
            output = {  'residueTag':int(item["begin"]),
                        'residueID': int(item["begin"]),
                        'residueName': resname,
                        'residueChain': "A",
                        'glycanIDs': get_filtered_glycanlist(config.data_dir,resname)
                            }           
            final.append(output)
            uniprot_spots.append(int(item["begin"]))
    outputSequence = uniprotSequence["sequence"]
    outputSequenceLength = uniprotSequence["length"]
    output = {
        "sequenceLength": outputSequenceLength,
        "sequence": outputSequence,
        "glycosylations": uniprotGlycosylations,
    }  
    sequence_with_info, sequences = pdb.get_sequence_from_pdb(outputFilePath)
    spots = pdb.find_glycosylation_spots(sequence_with_info)
    for i in spots:
        if i[1] not in uniprot_spots:
            resname = df.loc[(df['ResId'] == int(i[1])) & (df['Chain'] == i[2]), 'ResName'].iloc[0]
            output2 = {  
                        'residueTag': int(i[0]),
                        'residueID': int(i[1]),
                        'residueName': resname,
                        'residueChain': "A",
                        'glycanIDs':  get_filtered_glycanlist(config.data_dir,resname),
                            }  
            final.append(output2)
    return jsonify({'uniprot': uniprotID, 'glycosylation_locations': output, 'requestURL': requestURLcif, 'configuration' : final})


@app.route('/api/process_uniprot', methods=['POST'])
def process_uniprot():
    downloadLocation = config.upload_dir

    # Extract data from the JSON payload
    data = request.json
    uniprotID = data.get('uniprotID')
    selectedGlycans = data.get('selectedGlycans')

    # Extract glycans and glycosylation_locations from selectedGlycans
    glycans = [glycan for glycan in selectedGlycans.values() if glycan]  # Filter out empty glycans
    glycosylation_locations = [location for location in selectedGlycans.keys()]

    # Assuming you have the necessary pdb and algo modules and methods
    protein = pdb.parse(f'{downloadLocation}{uniprotID}.pdb')
    g, clash, box,link_pairs = algo.attach(protein, glycans, glycosylation_locations)
    
    # outputshortfilepath = f'{uniprotID}_glycosylated_{time.time()}.pdb'
    now = datetime.now()
    outputshortfilepath = f'{uniprotID.strip(".pdb")}_reglyco_{now.strftime("%Y%m%d%H%M")}.pdb'
    
    outputfilepath = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.exportPDB(outputfilepath, pdb.to_normal(g), link_pairs)
    
    return jsonify({'output':outputshortfilepath, 'clash': clash, 'box': box})


@app.route('/api/one_uniprot', methods=['POST'])
def one_uniprot():
    downloadLocation = config.upload_dir

    # Extract data from the JSON payload
    data = request.json
    uniprotID = data.get('uniprotID')
    outputFileName = uniprotID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
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

    selectedGlycans = {}

    for item in uniprotFeatures:
        if item["type"] == "CARBOHYD":
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
    g, clash, box,link_pairs = algo.attach(protein, glycans, glycosylation_locations)
    # outputshortfilepath = f'{uniprotID}_glycosylated_{time.time()}.pdb'
    now = datetime.now()
    outputshortfilepath = f'{uniprotID.strip(".pdb")}_reglyco_{now.strftime("%Y%m%d%H%M")}.pdb'
    
    outputfilepath = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.exportPDB(outputfilepath, pdb.to_normal(g),link_pairs)
    
    return jsonify({'output':outputshortfilepath, 'clash': clash, 'box': box})


@app.route('/api/upload_pdb', methods=['POST'])
def upload_pdb():
    upload_dir = config.upload_dir
    if 'pdbFile' not in request.files:
        return jsonify(error="No file part"), 400
    file = request.files['pdbFile']
    # Check if the post request has the file part
    if file.filename == '':
        return jsonify(error="No selected file"), 400

    if file :
        filename = os.path.join(upload_dir, file.filename)
        file.save(filename)
        pdb.remove_hydrogens(filename,filename)
        
        base_dir = config.data_dir
        if os.path.exists(base_dir):
            dir_list = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

        outputfilename = f'{filename}_glycosylated.pdb'
        protein = pdb.parse(filename)
        df = pdb.to_DF(protein)
        final = []
        sequence_with_info, sequences = pdb.get_sequence_from_pdb(filename)
        spots = pdb.find_glycosylation_spots(sequence_with_info)
        for i in spots:
             resname = df.loc[(df['ResId'] == int(i[1])) & (df['Chain'] == i[2]), 'ResName'].iloc[0]
             output = {  
                        'residueTag': int(i[0]),
                        'residueID': int(i[1]),
                        'residueName': resname,
                        'residueChain': i[2],
                        'glycanIDs':  get_filtered_glycanlist(config.data_dir,resname),
                            }  
             final.append(output)
        data = {
            "sequenceLength": len(sequence_with_info),
            "sequence": sequences ,
            "glycosylations": spots,
        }
        return jsonify({'uniprot':f'{file.filename}','glycosylation_locations': data, 'requestURL': f'https://glycoshape.io/output/{file.filename}', 'configuration' : final})
    return jsonify(error="Invalid file type"), 400

@app.route('/api/oneshot_pdb', methods=['POST'])
def oneshot_pdb():
    downloadLocation = config.upload_dir

    # Extract data from the JSON payload
    data = request.json
    customPDB = data.get('customPDB')
    if customPDB:
        filename = data.get('filename')
    else:
        filename = data.get('filename')+".pdb"
    result_residue = data.get('result')
    selectedGlycan = data.get('selectedGlycanOption')
     # Assuming you have the necessary pdb and algo modules and methods
    protein = pdb.parse(f'{downloadLocation}{filename}')
    selectedGlycans = {}
    for item in result_residue:
        if item.get('clash_solved', False):
            residue = item.get('residue', '')
            if residue:
                # Insert an underscore before the last character
                formatted_residue = residue[:-1] + '_' + residue[-1]
                selectedGlycans[formatted_residue] = selectedGlycan
    
    # Extract glycans and glycosylation_locations from selectedGlycans
    glycans = [glycan for glycan in selectedGlycans.values() if glycan]  # Filter out empty glycans
    glycosylation_locations = [location for location in selectedGlycans.keys()]
   
   
    g, clash , box ,link_pairs= algo.attach_skip(protein, glycans, glycosylation_locations)
    # outputshortfilepath = f'{filename}_glycosylated_{time.time()}.pdb'
    now = datetime.now()
    outputshortfilepath = f'{filename.strip(".pdb")}_reglyco_{now.strftime("%Y%m%d%H%M")}.pdb'
    
    outputfilepath = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.exportPDB(outputfilepath, pdb.to_normal(g),link_pairs)

    return jsonify({'output':outputshortfilepath, 'clash': clash, 'box': box})

@app.route('/api/process_pdb', methods=['POST'])
def process_pdb():
    downloadLocation = config.upload_dir

    # Extract data from the JSON payload
    data = request.json
    uniprotID = data.get('uniprotID')
    selectedGlycans = data.get('selectedGlycans')

    # Extract glycans and glycosylation_locations from selectedGlycans
    glycans = [glycan for glycan in selectedGlycans.values() if glycan]  # Filter out empty glycans
    glycosylation_locations = [location for location in selectedGlycans.keys()]
    # print(glycosylation_locations)
    # Assuming you have the necessary pdb and algo modules and methods
    protein = pdb.parse(f'{downloadLocation}{uniprotID}')
    g, clash , box,link_pairs = algo.attach(protein, glycans, glycosylation_locations)
    # outputshortfilepath = f'{uniprotID}_glycosylated_{time.time()}.pdb'
    now = datetime.now()
    outputshortfilepath = f'{uniprotID.strip(".pdb")}_reglyco_{now.strftime("%Y%m%d%H%M")}.pdb'
    
    
    outputfilepath = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.exportPDB(outputfilepath, pdb.to_normal(g),link_pairs)

    return jsonify({'output':outputshortfilepath, 'clash': clash, 'box': box})


@app.route('/api/scan', methods=['POST'])
def scan():
    downloadLocation = config.upload_dir

    # Extract data from the JSON payload
    data = request.json
    customPDB = data.get('customPDB')
    if customPDB:
        filename = data.get('filename')
    else:
        filename = data.get('filename')+".pdb"
    selectedGlycan = "GlcNAc"
     # Assuming you have the necessary pdb and algo modules and methods
    protein = pdb.parse(f'{downloadLocation}{filename}')
    df = pdb.to_DF(protein)
    sequence_with_info, sequences = pdb.get_sequence_from_pdb(f'{downloadLocation}{filename}')
    spots = pdb.find_glycosylation_spots_N(sequence_with_info)
    selectedGlycans = {}
    for i in spots:
            resname = df.loc[(df['ResId'] == int(i[1])) & (df['Chain'] == i[2]), 'ResName'].iloc[0]
            key = f"{int(i[1])}_{i[2]}" 
            selectedGlycans[key] = selectedGlycan
    
    # Extract glycans and glycosylation_locations from selectedGlycans
    glycans = [glycan for glycan in selectedGlycans.values() if glycan]  # Filter out empty glycans
    glycosylation_locations = [location for location in selectedGlycans.keys()]
   
    g, clash , box,link_pairs = algo.attach_skip(protein, glycans, glycosylation_locations)

    # Regex pattern to capture the relevant data
    pattern = re.compile(
        r"Residue :\s*(\d+[A-Z])\n\s*(.*?) has \d+ Clusters"
        r"(.*?)(?:Clash Solved for (.*?) at residue :\1 with phi :\s*(\d+) and psi :\s*(-?\d+), cluster (\d+)|Clash exist for (.*?) at residue :\1 with phi :\s*(\d+) and psi :\s*(-?\d+), cluster (\d+))",
        re.DOTALL
    )
    
    # Find all matches
    matches = pattern.findall(box)
    
    # Process matches to build the JSON structure
    results = []

    for match in matches:
        residue, glycan, _, solved_glycan, phi, psi, cluster, exist_glycan, exist_phi, exist_psi, exist_cluster = match
        result = {
            'residue': residue,
            'glycan': glycan or solved_glycan or exist_glycan,
            'clash_solved': bool(solved_glycan),
            'phi': int(phi or exist_phi),
            'psi': int(psi or exist_psi),
            'cluster': int(cluster or exist_cluster)
        }
        results.append(result)
    now = datetime.now()
    outputshortfilepath = f'{filename.strip(".pdb")}_reglyco_{now.strftime("%Y%m%d%H%M")}.pdb'
    outputfilepath = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.exportPDB(outputfilepath, pdb.to_normal(g),link_pairs)

    return jsonify({'output':outputshortfilepath, 'clash': clash, 'box': box, 'results': results})


@app.route('/api/maturation', methods=['POST'])
def maturation():
    downloadLocation = config.upload_dir

    # Extract data from the JSON payload
    data = request.json
    customPDB = data.get('customPDB')
    if customPDB:
        filename = data.get('filename')
    else:
        filename = data.get('filename')+".pdb"
    residue = data.get('selectedResidueMaturation')
    selectedGlycan = "GlcNAc"

    
     # Assuming you have the necessary pdb and algo modules and methods
    protein = pdb.parse(f'{downloadLocation}{filename}')
    df = pdb.to_DF(protein)
    sequence_with_info, sequences = pdb.get_sequence_from_pdb(f'{downloadLocation}{filename}')
    spots = pdb.find_glycosylation_spots_N(sequence_with_info)
    selectedGlycans = {}
    for i in spots:
            resname = df.loc[(df['ResId'] == int(i[1])) & (df['Chain'] == i[2]), 'ResName'].iloc[0]
            key = f"{int(i[1])}_{i[2]}" 
            selectedGlycans[key] = selectedGlycan
    
    # Extract glycans and glycosylation_locations from selectedGlycans
    glycans = [glycan for glycan in selectedGlycans.values() if glycan]  # Filter out empty glycans
    glycosylation_locations = [location for location in selectedGlycans.keys()]
   
   
    g, clash , box , link_pairs = algo.attach_skip(protein, glycans, glycosylation_locations)

    # Regex pattern to capture the relevant data
    pattern = re.compile(
        r"Residue :\s*(\d+[A-Z])\n\s*(.*?) has \d+ Clusters"
        r"(.*?)(?:Clash Solved for (.*?) at residue :\1 with phi :\s*(\d+) and psi :\s*(-?\d+), cluster (\d+)|Clash exist for (.*?) at residue :\1 with phi :\s*(\d+) and psi :\s*(-?\d+), cluster (\d+))",
        re.DOTALL
    )
    
    # Find all matches
    matches = pattern.findall(box)
    
    # Process matches to build the JSON structure
    results = []

    for match in matches:
        residue, glycan, _, solved_glycan, phi, psi, cluster, exist_glycan, exist_phi, exist_psi, exist_cluster = match
        result = {
            'residue': residue,
            'glycan': glycan or solved_glycan or exist_glycan,
            'clash_solved': bool(solved_glycan),
            'phi': int(phi or exist_phi),
            'psi': int(psi or exist_psi),
            'cluster': int(cluster or exist_cluster)
        }
        results.append(result)
    
    now = datetime.now()
    outputshortfilepath = f'{filename.strip(".pdb")}_reglyco_{now.strftime("%Y%m%d%H%M")}.pdb'
    outputfilepath = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.exportPDB(outputfilepath, pdb.to_normal(g),link_pairs)

    return jsonify({'output':outputshortfilepath, 'clash': clash, 'box': box, 'results': results})

@app.route('/api/log', methods=['GET'])
def log_visitor():
    # Check for existing cookie
    if request.cookies.get('visited'):
        return 'Already logged today', 200

    # Get client IP address
    if request.headers.getlist("X-Forwarded-For"):
        ip = request.headers.getlist("X-Forwarded-For")[0]
    else:
        ip = request.remote_addr

    # Log IP and timestamp
    timestamp = datetime.now()
    with open("visitors.log", "a") as log_file:
        log_file.write(f"{timestamp}: {ip}\n")

    # Create a response and set a cookie with 24-hour expiration
    response = make_response("Logged", 200)
    response.set_cookie('visited', 'yes', max_age=86400)  # 86400 seconds = 24 hours
    return response



@app.route('/api/upload_pdb_swap', methods=['POST'])
def upload_pdb_swap():
    upload_dir = config.upload_dir
    if 'pdbFile' not in request.files:
        return jsonify(error="No file part"), 400
    file = request.files['pdbFile']
    # Check if the post request has the file part
    if file.filename == '':
        return jsonify(error="No selected file"), 400

    if file :
        filename = os.path.join(upload_dir, file.filename)
        file.save(filename)
        pdb.remove_hydrogens(filename,filename)
        
        base_dir = config.data_dir
        if os.path.exists(base_dir):
            dir_list = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

        outputfilename = f'{filename}_glycosylated.pdb'
        protein = pdb.parse(filename)
        df = pdb.to_DF(protein)
        final = []
        sequence_with_info, sequences = pdb.get_sequence_from_pdb(filename)
        spots = pdb.find_glycosylation_spots_N(sequence_with_info)
        for i in spots:
             resname = df.loc[(df['ResId'] == int(i[1])) & (df['Chain'] == i[2]), 'ResName'].iloc[0]
             output = {  
                        'residueTag': int(i[0]),
                        'residueID': int(i[1]),
                        'residueName': resname,
                        'residueChain': i[2],
                        'glycanIDs':  [True,False],
                            }  
             final.append(output)
        data = {
            "sequenceLength": len(sequence_with_info),
            "sequence": sequences ,
            "glycosylations": spots,
        }
        return jsonify({'uniprot':f'{file.filename}','glycosylation_locations': data, 'requestURL': f'https://glycoshape.io/output/{file.filename}', 'configuration' : final})
    return jsonify(error="Invalid file type"), 400


@app.route('/api/process_pdb_swap', methods=['POST'])
def process_pdb_swap():
    downloadLocation = config.upload_dir

    # Extract data from the JSON payload
    data = request.json
    uniprotID = data.get('uniprotID')
    selectedGlycans = data.get('selectedGlycans')
    selected_list = [location for location in selectedGlycans.keys()]
    

    now = datetime.now()
    outputshortfilepath = f'{uniprotID.strip(".pdb")}_swapped_{now.strftime("%Y%m%d%H%M")}.pdb'
    outputfilepath = f'{downloadLocation}{outputshortfilepath}'
    
    pdb.swap_residues(f'{downloadLocation}{uniprotID}',selected_list,outputfilepath)

    return jsonify({'output':outputshortfilepath, 'clash': False, 'box': "No error."})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000)
