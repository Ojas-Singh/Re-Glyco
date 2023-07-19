# Information

N-linked
Phi = CG-ND2-C1-O5
Psi = CB-CG-ND2-C1

O-linked
Phi = CB-OG1-C1-O5
Psi = CA-CB-OG1-C1

C-linked
Phi = CG-CD1-C1-O5


# Installation


```


conda create -n reglyco python=3.10
conda activate reglyco
pip install -r requirements.txt

sudo apt install build-essential

cd glycors
maturin develop --release
cd ..
streamlit run main.py
```

python3 -m pip -r requirements.txt
export PATH="$HOME/.local/bin:$PATH"


modify config.py 
run

python main.py

gunicorn -w 4 api.py:app

flask -A api.py run     


# things for oracle for Hosting Websites
```
sudo iptables -P INPUT ACCEPT
sudo iptables -P OUTPUT ACCEPT
sudo iptables -P FORWARD ACCEPT
sudo iptables -F
sudo ufw allow 22
sudo ufw allow 80
sudo ufw allow 8080
sudo ufw allow 443
sudo ufw enable
```

localhost/loopback
```
sudo iptables -t nat -I OUTPUT -p tcp -d 127.0.0.1 --dport 80 -j REDIRECT --to-ports 3000
```
external
```
sudo iptables -t nat -I PREROUTING -p tcp --dport 80 -j REDIRECT --to-ports 3000
```
```
screen -S name
screen -r name
ctr a + d   -> to detach
pkill screen

```


# Re-Glyco API Documentation

This document provides the details needed to utilize the API endpoints for GlycoShape, a service for protein glycosylation modeling. Please note that this API is currently unreleased and subject to change. The server address is `glycoshape.healoor.me:8000`.

## /available_glycans - GET

This endpoint retrieves all available glycans. No parameters are required for this GET request.

Example:
```
curl -X GET glycoshape.healoor.me:8000/available_glycans
```

## /linked_glycans - POST

This endpoint retrieves a list of glycans based on the link type. It expects a JSON body with one parameter `link_type`, which can be 'N_linked', 'O_linked', 'C_linked', or 'X_linked'.

Example:
```bash
curl -X POST -H "Content-Type: application/json" -d '{"link_type":"N_linked"}' glycoshape.healoor.me:8000/linked_glycans
```

## /custom_pdb_spots - POST

This endpoint accepts a PDB file and returns potential glycosylation spots. The request should include a file with key 'file'.

Example:
```bash
curl -X POST -F "file=@/path/to/your/file.pdb" glycoshape.healoor.me:8000/custom_pdb_spots
```

Replace `@/path/to/your/file.pdb` with the path to your file.

## /process_pdb - POST

This endpoint accepts a PDB file and a list of glycans and glycosylation locations, then returns a processed PDB file (in base64 format) and clash information. The request should include a file with key 'file', and the 'glycans' and 'glycosylation_locations' data should be provided as form data, not JSON. 

Example:
```bash
curl -X POST -F "file=@/path/to/your/file.pdb" -F "glycans=DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH,DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH" -F "glycosylation_locations=38,75" glycoshape.healoor.me:8000/process_pdb
```

Replace `@/path/to/your/file.pdb` with the path to your file. The `glycans` and `glycosylation_locations` data should be comma-separated lists of strings and integers, respectively.

The returned `file` data will be a base64 encoded string. To decode this and write it to a file in Python, you can do:

```python
import base64

data = "base64 encoded string"  # replace this with the actual string
decoded_data = base64.b64decode(data)
with open("output.pdb", "wb") as f:
    f.write(decoded_data)
```

```python
import requests
import json
import base64

data = {
    'glycans': ['DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH','DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH'],  # Replace with the actual glycans
    'glycosylation_locations': [38,75]  # Replace with the actual glycosylation locations
}
# Open the PDB file in binary mode
with open('output/FOLD.pdb', 'rb') as f:
    file = {'file': f}

    # Send the request
    response = requests.post('http://glycoshape.healoor.me:8000/process_pdb', files=file, data=data)

# Parse the response JSON
response_data = response.json()

# The 'file' field in the response contains the base64-encoded PDB file
encoded_file = response_data['file']

# The 'clash' field in the response is the clash info
clash = response_data['clash']

# Decode the base64-encoded PDB file
decoded_file = base64.b64decode(encoded_file)

# Write the decoded PDB file to disk
with open('output.pdb', 'wb') as f:
    f.write(decoded_file)

print(f"Clash info: {clash}")

```

## /process_uniprot - POST

This endpoint accepts a Uniprot ID, a list of glycans, and glycosylation locations. It returns a processed PDB file (in base64 format) and clash information. The request should include a JSON body with keys 'uniprot', 'glycans', and 'glycosylation_locations'.

Example:
```bash
curl -X POST -H "Content-Type: application/json" -d '{"uniprot":"P29016","glycans":["DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH","DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH"],"glycosylation_locations":[38,75]}' glycoshape.healoor.me:8000/process_uniprot
```

Again, the returned `file` data will be a base64 encoded string, which can be decoded and written to a file as shown in the previous section.

```python
import requests
import json
import base64

# Define the URL of your Flask application
url = "http://glycoshape.healoor.me:8000/process_uniprot"

# Define the parameters to be sent in the request
data = {
    'uniprot': 'P29016',  # Replace with the actual uniprot ID
    'glycans': ['DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH','DManpb1-4DGlcpNAcb1-4DGlcpNAca1-OH'],  # Replace with the actual glycans
    'glycosylation_locations': [38,75]  # Replace with the actual glycosylation locations
}

# Make the POST request
response = requests.post(url, json=data)

# Parse the response JSON
response_data = response.json()

# The 'file' field in the response contains the base64-encoded PDB file
encoded_file = response_data['file']

# The 'clash' field in the response is the clash info
clash = response_data['clash']

# Decode the base64-encoded PDB file
decoded_file = base64.b64decode(encoded_file)

# Write the decoded PDB file to disk
with open('output.pdb', 'wb') as f:
    f.write(decoded_file)

print(f"Clash info: {clash}")
```

