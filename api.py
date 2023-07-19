from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from tempfile import NamedTemporaryFile
import os
import time
import config
import base64
from lib import pdb,algo,prep

app = Flask(__name__)
CORS(app)

def get_glycan_list(suffix):
    directory = config.data_dir
    folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder)) and folder.endswith(suffix)]
    return folders

@app.route('/available_glycans', methods=['GET'])
def list_glycans():
    base_dir = config.data_dir
    if os.path.exists(base_dir):
        dir_list = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
        return jsonify({'glycan_list': dir_list})
    else:
        return jsonify({'error': 'Directory not found'})
    
@app.route('/linked_glycans', methods=['POST'])
def linked_glycans():
    link_type = request.json.get('link_type')
    if link_type == 'N_linked':
        dir_list = ["None",] + get_glycan_list("DGlcpNAca1-OH") + get_glycan_list("DGlcpNAcb1-OH")
    elif link_type == 'O_linked':
        dir_list = ["None",] + get_glycan_list("DGalpNAca1-OH") + get_glycan_list("DGalpNAcb1-OH")
    elif link_type == 'C_linked':
        dir_list = ["None",] + get_glycan_list("DManpa1-OH") + get_glycan_list("DManpb1-OH")
    elif link_type == 'X_linked':
        dir_list = ["None",] + get_glycan_list("DXylpa1-OH") + get_glycan_list("DXylpb1-OH")
    else:
        return jsonify({'error': 'Invalid link type'})

    return jsonify({'glycan_list': dir_list})

@app.route('/custom_pdb_spots', methods=['POST'])
def custom_pdb_spots():
    file = request.files['file']
    if not file:
        return "No file", 400

    temp_file = NamedTemporaryFile(delete=False)
    file.save(temp_file.name)
    temp_file.close()
    protein = pdb.parse(temp_file.name)
    protein_df= pdb.to_DF(protein)
    spots= protein_df.loc[(protein_df['ResName']=="ASN") & (protein_df['Name']== 'CB') |(protein_df['ResName']=="THR") | (protein_df['ResName']=="TYR")|((protein_df['ResName']=="TRP") | (protein_df['ResName']=="SER")) & (protein_df['Name']== 'CB') ,['ResId']].iloc[:]['ResId'].tolist()
    os.unlink(temp_file.name)

    return jsonify({'spots': spots})

@app.route('/process_pdb', methods=['POST'])
def process_pdb():
    file = request.files['file']

    # Use request.form.get instead of request.json.get
    glycans = request.form.get('glycans').split(',')
    glycosylation_locations = [int(x) for x in request.form.get('glycosylation_locations').split(',')]

    if not file:
        return "No file", 400

    with NamedTemporaryFile(delete=False) as temp_file:
        file.save(temp_file.name)
        temp_file_path = temp_file.name  # Store the file path to delete it later

    protein = pdb.parse(temp_file_path)
    g, clash = algo.attach(protein, glycans, glycosylation_locations)
    os.unlink(temp_file_path)  # Delete the original temp file

    with NamedTemporaryFile(suffix=".pdb", delete=False) as temp_file:
        pdb.exportPDB(temp_file.name, pdb.to_normal(g))

        # The file will be closed after the with block, no need for f.close()
        with open(temp_file.name, "rb") as f:
            encoded_string = base64.b64encode(f.read()).decode('utf-8')

    os.unlink(temp_file.name)  # Don't forget to delete the temp file after using it

    return jsonify({'file': encoded_string, 'clash': clash})

@app.route('/process_uniprot', methods=['POST'])
def process_uniprot():
    uniprot = request.json.get('uniprot')
    glycans = request.json.get('glycans')
    glycosylation_locations = request.json.get('glycosylation_locations')
    file_path = prep.download_and_prepare_alphafoldDB_model(uniprot, "output/temp/")
    protein = pdb.parse(file_path)
    g, clash = algo.attach(protein, glycans, glycosylation_locations)

    temp_file = NamedTemporaryFile(suffix=".pdb", delete=False)
    pdb.exportPDB(temp_file.name, pdb.to_normal(g))

    with open(temp_file.name, "rb") as f:
        encoded_string = base64.b64encode(f.read()).decode('utf-8')
    f.close()
    # os.unlink(temp_file.name)

    return jsonify({'file': encoded_string, 'clash': clash})



if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)
