data_dir = "DB/"
upload_dir = "output/"
domain_name = "https://glycoshape.io"
N_linked = {"Res":["ASN"],"phi": (-130,-63),"psi":(152,205)}
O_linked = {"Res":["SER","THR"],"phi": (65,86),"psi":(111,164)}
C_linked = {"Res":["TYR","TRP"],"phi": (110,150),"psi":(-3,3)}



ASN = {"link": "beta"  , "A": ["CB"], "B": ["CG"], "C": ["ND2"], "D": ["C1"], "E": ["O5"],
        "sugars": {
        "GlcNAc": {"phi": (-130,-63), "psi":(152,205),"link": "beta"},
        }
}
THR = {"link": "alpha" , "A": ["CA"], "B": ["CB"], "C": ["OG1","OG"], "D": ["C1"], "E": ["O5"],
        "sugars": {
        "GalNAc": {"phi": (55,83), "psi":(86,142),"link": "alpha"},
        "Fuc": {"phi": (274,304), "psi":(148,182),"link": "alpha"},
        "Man": {"phi": (61,88), "psi":(88,150),"link": "alpha"},
        "GlcNAc": {"phi": (143,221), "psi":(171,193),"link": "alpha"},
        }
}
       
SER = {"link": "alpha" , "A": ["CA"], "B": ["CB"], "C": ["OG1","OG"], "D": ["C1"], "E": ["O5"],
        "sugars": {
        "GalNAc": {"phi": (61,86), "psi":(163,224),"link": "alpha"},
        "Fuc": {"phi": (58,106), "psi":(86,192),"link": "alpha"},
        "Glc": {"phi": (261,297), "psi":(146,219),"link": "beta"},
        "Xyl": {"phi": (265,309), "psi":(107,231),"link": "beta"},
        "GlcNAc": {"phi": (272,308), "psi":(178,300),"link": "beta"},
        }
}
       
TYR = {"link": "alpha" , "A": ["CB"], "B": ["CG"], "C": ["CD1"], "D": ["C1"], "E": ["O5"],
        "sugars": {
        "Man": {"phi": (110,150), "psi":(-3,-3),"link": "alpha"}},
        
}
       
TRP = {"link": "alpha" , "A": ["CB"], "B": ["CG"], "C": ["CD1"], "D": ["C1"], "E": ["O5"],
        "sugars": {
        "Man": {"phi": (110,150), "psi":(-3,3),"link": "alpha"},
        }
}

RESIDUE_MAP = {
    "ASN": ASN,
    "THR": THR,
    "SER": SER,
    "TRP": TRP
}