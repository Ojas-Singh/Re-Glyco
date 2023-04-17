from stmol import showmol
from lib import pdb
import py3Dmol
import streamlit as st

@st.cache_data
def show3d(fold,system,confidence,glycosylation_locations):
    protein = pdb.parse(fold)
    xyzview = py3Dmol.view()
    xyzview.addModelsAsFrames(system)
    colors = ["#FF7D45" for x in range(50)]+["#FFDB13" for x in range(20)]+["#65CBF3" for x in range(20)]+["#0053D6" for x in range(15)]
    xyzview.setStyle({'cartoon':{'color':'grey'}})
    for resid in range(len(confidence)):
        xyzview.addStyle({"chain": "A", "resi": str(resid+1)},
                        {"cartoon": {"color":colors[int(confidence[resid])]}})
    xyzview.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "lightgrey"},{"hetflag": False})
    for i in range(len(glycosylation_locations)):
        xyzview.addStyle({"chain": "A", "resi": glycosylation_locations[i]["begin"], "elem": "C"},
                        {"stick": {"color": colors[int(confidence[int(glycosylation_locations[i]["begin"])-1])], "radius":  0.2}})
        xyzview.addStyle({"chain": "A", "resi": glycosylation_locations[i]["begin"]},
                            {"stick": {"radius":  0.4}})    
        xyzview.addResLabels({"chain": "A","resi": glycosylation_locations[i]["begin"]},
        {"backgroundColor": "lightgray","fontColor": "purple","backgroundOpacity": 0.5})
    xyzview.setBackgroundColor('#FFFFFF')
    xyzview.zoomTo()
    showmol(xyzview,height=400,width=800)


def show3doutput(g1,glycosylation_locations,confidence):
    colors = ["#FF7D45" for x in range(50)]+["#FFDB13" for x in range(20)]+["#65CBF3" for x in range(20)]+["#0053D6" for x in range(15)]
    xyzview1 = py3Dmol.view()
    xyzview1.addModelsAsFrames(g1)
    for n,chain,color in zip(range(len(glycosylation_locations)+1),list("ABCDEFGH"),
                                ["grey","#FF7B89","#8A5082","#6F5F90","#758EB7","#A5CAD2","blue","orange"]):
                    if chain=="A":
                        xyzview1.setStyle({'chain':chain},{'stick': {"opacity": 0.6,'color':color,"radius":  0.2}})
                        xyzview1.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "lightgrey"},{"hetflag": False})
                    else:
                        xyzview1.setStyle({'chain':chain},{'stick': {'color':color, "radius":  0.4}})
    for resid in range(len(confidence)):
        xyzview1.addStyle({"chain": "A", "resi": str(resid+1)},
                        {"stick": {"color":colors[int(confidence[resid])], "radius":  0.4}})
    for i in range(len(glycosylation_locations)):
        xyzview1.addStyle({"chain": "A", "resi": glycosylation_locations[i]["begin"], "elem": "C"},
                        {"stick": {"color": "red", "radius":  0.2}})

        xyzview1.addStyle({"chain": "A", "resi": glycosylation_locations[i]["begin"]},
                            {"stick": {"radius":  0.4}})
    xyzview1.setBackgroundColor('#FFFFFF')
    xyzview1.zoomTo()
    showmol(xyzview1,height=400,width=1100)