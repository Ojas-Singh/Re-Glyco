#demo code for uploading pdb file and downloading processed pdb file

import requests
import os
import json
from datetime import datetime

# Constants
API_BASE_URL = "https://glycoshape.org"  # API URL
# API_BASE_URL = "http://127.0.0.1:8000" # API URL
UPLOAD_ENDPOINT = "/api/upload_pdb"
PROCESS_ENDPOINT = "/api/process_pdb"
UPLOAD_DIR = "output"  # Local directory to save downloaded PDB files

# Ensure upload directory exists
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR)

def upload_pdb(file_path):
    """Uploads a PDB file and gets glycosylation configurations."""
    with open(file_path, 'rb') as file:
        files = {'pdbFile': file}
        response = requests.post(API_BASE_URL + UPLOAD_ENDPOINT, files=files)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Error uploading PDB file: {response.json()}")

def process_pdb(uniprot_id, glycan_configurations):
    """Processes a PDB file with given glycan configurations."""
    data = {
        'uniprotID': uniprot_id,
        'selectedGlycans': glycan_configurations
    }
    try:
        # Using a more realistic timeout value (e.g., 5 minutes)
        response = requests.post(API_BASE_URL + PROCESS_ENDPOINT, json=data, timeout=4000)
        response.raise_for_status()  # Raises an HTTPError for bad responses
        return response.json()
    except requests.exceptions.Timeout:
        # Handle request timeout separately
        raise TimeoutError("The request timed out. Please try again later.")
    except requests.exceptions.HTTPError as err:
        # Handle HTTP errors
        error_msg = response.text  # Fallback to response text if JSON is not available
        try:
            error_json = response.json()
            error_msg = error_json.get('message', error_msg)  # Attempt to get a more detailed message
        except ValueError:
            pass  # Use the default error_msg if JSON parsing fails
        raise Exception(f"Error processing PDB file: {error_msg}") from err
    except requests.exceptions.RequestException as err:
        # Handle other requests exceptions
        raise Exception("An error occurred while processing the PDB file.") from err
def download_processed_pdb(file_name):
    """Downloads the processed PDB file."""
    response = requests.get(f"{API_BASE_URL}/output/{file_name}")
    if response.status_code == 200:
        with open(os.path.join(UPLOAD_DIR, file_name), 'wb') as file:
            file.write(response.content)
    else:
        raise Exception("Error downloading processed PDB file")

def main():
    # Example usage
    pdb_file_path = 'AF-P29016-F1-model_v4.pdb'  # Replace with actual PDB file path
    upload_response = upload_pdb(pdb_file_path)


    glycan_configurations = {
        '38_A': "Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",   # ---> residueID_residueChain :  glycanID (glycanID of choice from corresponding configurations)
        '75_A': "Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
        '146_A': "Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    }


    file_id = upload_response['uniprot']
    print(f"Uploaded PDB file: {pdb_file_path}")
    process_response = process_pdb(file_id, glycan_configurations)
    output_file_name = process_response['output']
    try:
        download_processed_pdb(output_file_name)
    except Exception as e:
        print(f"Error downloading processed PDB file moslty running locally check temp_files for output strucutres")
    print(f"Processed PDB file downloaded: {output_file_name}")
    print(process_response['box'])


if __name__ == "__main__":
    main()