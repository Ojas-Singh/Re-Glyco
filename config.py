data_dir = "/mnt/database/DB/"
upload_dir = "/mnt/database/server_dir/"
domain_name = "http://glycoshape.healoor.me:8088/"
N_linked = {"Res":["ASN"],"phi": (-130,-63),"psi":(152,205)}
O_linked = {"Res":["SER","THR"],"phi": (65,86),"psi":(111,164)}
C_linked = {"Res":["TYR","TRP"],"phi": (110,150),"psi":(-3,3)}



ASN = {"link": "beta"  , "A": ["CB"], "B": ["CG"], "C": ["ND2"], "D": ["C1"], "E": ["O5"], "phi": (-130,-63), "psi":(152,205)}
THR = {"link": "alpha" , "A": ["CA"], "B": ["CB"], "C": ["OG1","OG"], "D": ["C1"], "E": ["O5"], "phi": (65,86), "psi":(111,164)}
SER = {"link": "alpha" , "A": ["CA"], "B": ["CB"], "C": ["OG1","OG"], "D": ["C1"], "E": ["O5"], "phi": (65,86), "psi":(111,164)}
TYR = {"link": "alpha" , "A": ["CB"], "B": ["CG"], "C": ["CD1"], "D": ["C1"], "E": ["O5"], "phi": (110,150), "psi":(-3,3)}
TRP = {"link": "alpha" , "A": ["CB"], "B": ["CG"], "C": ["CD1"], "D": ["C1"], "E": ["O5"], "phi": (110,150), "psi":(-3,3)}

RESIDUE_MAP = {
    "ASN": ASN,
    "THR": THR,
    "SER": SER,
    "TYR": TYR,
    "TRP": TRP
}