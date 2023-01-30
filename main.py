import streamlit as st
import pandas as pd
import numpy as np
from lib import prep,algo,pdb
import py3Dmol
from stmol import showmol
from PIL import Image
import time
from io import StringIO
import io
from annotated_text import annotated_text
import os

st.set_page_config(page_title="Re-GLyco", page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 
image = Image.open('logo.png')
st.sidebar.image(image, caption='')
uniprot=""
uni_id = st.sidebar.text_input(
        "Enter UNIPROT ID of protein...",placeholder ="P29016 etc"
    )

st.sidebar.write("Or upload PDB manually.")
uploaded_file = st.sidebar.file_uploader("")

@st.cache
def fetch(uni_id):
    fold = prep.download_and_prepare_alphafoldDB_model(uni_id,"output/temp/")
    out= prep.query_uniprot_for_glycosylation_locations(uni_id)
    return fold,out

if not uni_id=="" and uploaded_file is None:
    fold, out = fetch(uni_id)
    with open(fold) as ifile:
        system = "".join([x for x in ifile])
        confidence= []
        p=1
        lines = system.split("\n")
        for x in lines:
            if x.startswith("ATOM"):
                if int((x[23:27]).strip(" "))==p:
                    confidence.append(float((x[61:67]).strip(" ")))
                    p+=1
        tab1, tab2 = st.tabs(["Structure", "glycosylations"])
        glycosylation_locations = out["glycosylations"]
        with tab1:
            col1, col2 = st.columns(2)
            with col1:
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
            with col2:
                btn = st.download_button(
                        label="Download AlphaFold Structure",
                        data=system,
                        file_name=uni_id+".pdb",
                        mime='text/csv'
                    )
                container = st.container()
                container.write("Model Confidence:")
                annotated_text(("Very high", "(pLDDT > 90)", "#0053D6"))
                annotated_text(("Confident", "(90 > pLDDT > 70)", "#65CBF3"))
                annotated_text(("Low", "(70 > pLDDT > 50)", "#FFDB13"))
                annotated_text( ("Very low", "(pLDDT < 50)", "#FF7D45"))
                st.write("AlphaFold produces a per-residue confidence score (pLDDT) between 0 and 100. Some regions with low pLDDT may be unstructured in isolation. ")
        with tab2:
            for i in glycosylation_locations:
                st.write(i)
        glycans=[1 for x in range(len(glycosylation_locations))]
        glycosylation_locations_N=[]
        for i in range(len(glycosylation_locations)):
            if not glycosylation_locations[i]["description"].startswith('N-linked'):
                st.write("Only N Glycosylation yet! Spot :",glycosylation_locations[i]["begin"]," is ",glycosylation_locations[i]["description"])
            else:
                glycosylation_locations_N.append(glycosylation_locations[i])
        root='data'
        dirlist = [ item for item in os.listdir(root) if os.path.isdir(os.path.join(root, item)) ]
        for i in range(len(glycosylation_locations_N)): 
            options = st.selectbox(
                        'Which Glycans to attach on spot : '+str(glycosylation_locations_N[i]["begin"])+" ?" ,
                        (dirlist),key=str(i))
            glycans[i]=options
        if st.button('Process',key="process"):
            st.write('')
            s=time.time()
            with st.spinner('Processing...'):
                g = algo.attachwithwiggle(protein,glycans,glycosylation_locations_N)
            st.balloons()
            st.success("exec time : "+ str(int(time.time()-s)) +" seconds")
            g1 = pdb.exportPDB('output/out.pdb',pdb.to_normal(g))
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
            with open('output/out.pdb') as ofile:
                system = "".join([x for x in ofile])
                btn = st.download_button(
                    label="Download Re-Glyco Structure",
                    data=system,
                    file_name=uni_id+"_glycosylated.pdb",
                    mime='text/csv'
                )
        else:
            st.write('')
elif uploaded_file is not None:
    bytes_data = uploaded_file.getvalue()
    newFile = open("output/temp/custom.pdb", "wb")
    newFile.write(bytes_data)
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    string_data = stringio.read()
    # st.write(string_data)
    confidence= []
    p=1
    lines = string_data.split("\n")
    for x in lines:
        if x.startswith("ATOM"):
            if int((x[23:27]).strip(" "))==p:
                confidence.append(float((x[61:67]).strip(" ")))
                p+=1
    protein = pdb.parse("output/temp/custom.pdb")
    col1, col2 = st.columns(2)
    with col1:
        xyzview = py3Dmol.view()
        xyzview.addModelsAsFrames(string_data)
        colors = ["#FF7D45" for x in range(50)]+["#FFDB13" for x in range(20)]+["#65CBF3" for x in range(20)]+["#0053D6" for x in range(15)]
        xyzview.setStyle({'cartoon':{'color':'grey'}})
        for resid in range(len(confidence)):
            xyzview.addStyle({"chain": "A", "resi": str(resid+1)},
                            {"cartoon": {"color":colors[int(confidence[resid])]}})
        xyzview.setBackgroundColor('#FFFFFF')
        xyzview.zoomTo()
        showmol(xyzview,height=400,width=2000)
    with col2:
        container = st.container()
        container.write("Model Confidence:")
        annotated_text(("Very high", "(pLDDT > 90)", "#0053D6"))
        annotated_text(("Confident", "(90 > pLDDT > 70)", "#65CBF3"))
        annotated_text(("Low", "(70 > pLDDT > 50)", "#FFDB13"))
        annotated_text( ("Very low", "(pLDDT < 50)", "#FF7D45"))
    protein_df= pdb.to_DF(protein)
    spots= protein_df.loc[(protein_df['ResName']=="ASN") & (protein_df['Name']== 'CB'),['ResId']].iloc[:]['ResId'].tolist()
    # st.write()
    glycosylation_locations = st.multiselect(
    'Select Glycosylation Locations',
    spots,
    )
    glycans=[1 for x in range(len(glycosylation_locations))]
    for i in range(len(glycosylation_locations)):
        options = st.selectbox(
                    'Which Glycans to attach on spot : '+str(glycosylation_locations[i])+" ?" ,
                    ('bisecting','A3G3S1-F3', 'a2',"a2g2","a3g3","m5","m6_1","m6_2","m6_3","m7_1","m7_2","m7_3","m7_4","m8_1","m8_2","m8_3","m9"),key=str(i))
        glycans[i]=options
    if st.button('Process',key="process"):
            st.write('')
            s=time.time()
            with st.spinner('Processing...'):
                g = algo.attach(protein,glycans,glycosylation_locations)
            st.balloons()
            st.success("exec time :"+str(time.time()-s)+"Seconds")
            g1 = pdb.exportPDB('output/out.pdb',pdb.to_normal(g))
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
                xyzview1.addStyle({"chain": "A", "resi": glycosylation_locations[i], "elem": "C"},
                                {"stick": {"color": "red", "radius":  0.2}})
                xyzview1.addStyle({"chain": "A", "resi": glycosylation_locations[i]},
                                    {"stick": {"radius":  0.4}})
            xyzview1.setBackgroundColor('#FFFFFF')
            xyzview1.zoomTo()
            showmol(xyzview1,height=400,width=1100)
            with open('output/out.pdb') as ofile:
                system = "".join([x for x in ofile])
                btn = st.download_button(
                    label="Download Re-Glyco Structure",
                    data=system,
                    file_name=uni_id+"_glycosylated.pdb",
                    mime='text/csv'
                )
else:
    st.title('Just Another Glycosylator for Nerds')
    st.header('Under Construction!')
    st.subheader("")
    st.markdown("O15552")
    st.markdown("P29016")
    st.markdown("Q9BXJ4")

