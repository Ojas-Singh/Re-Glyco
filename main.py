import streamlit as st
import pandas as pd
import numpy as np
from lib import prep,algo,pdb
import py3Dmol
from stmol import showmol
from PIL import Image
import time

st.set_page_config(page_title="Re-GLyco", page_icon=None, layout="centered", initial_sidebar_state="auto", menu_items=None)

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
        "Enter UNIPROT ID of protien...",placeholder ="O15552 etc"
    )

st.sidebar.write("Or upload PDB manually.")
uploaded_file = st.sidebar.file_uploader("")

@st.cache
def fetch(uni_id):
    fold = prep.download_and_prepare_alphafoldDB_model(uni_id,"output/temp/")
    out= prep.query_uniprot_for_glycosylation_locations(uni_id)
    return fold,out


if not uni_id=="":
    fold, out = fetch(uni_id)
    with open(fold) as ifile:
        system = "".join([x for x in ifile])
        tab1, tab2, tab3 = st.tabs(["Structure", "glycosylations", "Owl"])
        glycosylation_locations = out["glycosylations"]
        with tab1:
            col1, col2 = st.columns(2)
            with col1:
                pass
            protein = pdb.parse(fold)
            xyzview = py3Dmol.view()
            xyzview.addModelsAsFrames(system)
            xyzview.setStyle({'cartoon':{'color':'grey'}})
            xyzview.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "lightgrey"},{"hetflag": False})

            for i in range(len(glycosylation_locations)):
                xyzview.addStyle({"chain": "A", "resi": glycosylation_locations[i]["begin"], "elem": "C"},
                                {"stick": {"color": "red", "radius":  0.2}})

                xyzview.addStyle({"chain": "A", "resi": glycosylation_locations[i]["begin"]},
                                    {"stick": {"radius":  0.4}})
                    
                xyzview.addResLabels({"chain": "A","resi": glycosylation_locations[i]["begin"]},
                {"backgroundColor": "lightgray","fontColor": "purple","backgroundOpacity": 0.5})

            xyzview.setBackgroundColor('#FFFFFF')
            xyzview.zoomTo()
            showmol(xyzview,height=400,width=500)
            with col2:
                btn = st.download_button(
                        label="Download AlphaFold Structure",
                        data=system,
                        file_name=uni_id+".pdb",
                        mime='text/csv'
                    )
        
        with tab2:
            for i in glycosylation_locations:
                st.write(i)
        glycans=[1 for x in range(len(glycosylation_locations))]
        for i in range(len(glycosylation_locations)):
            options = st.selectbox(
                        'What Glycans to attach? on spot '+str(glycosylation_locations[i]["begin"]) ,
                        ('bisecting','A3G3S1-F3', 'a2',"a2g2","a3g3","m5","m6_1","m6_2","m6_3","m7_1","m7_2","m7_3","m7_4","m8_1","m8_2","m8_3","m9"),key=str(i))
            # 
            if options=='A3G3S1-F3':
                picture = Image.open('data/'+options+'.png')
                st.image(picture, caption='Neu5Ac(a3)Gal(b4)[Fuc(a3)]GlcNAc(b6)[Gal(b4)GlcNAc(b2)]Man(a6)[Gal(b4)GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)GlcNAc()',width=300)
            glycans[i]=options
        print(glycans)
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
                                xyzview1.setStyle({'chain':chain},{'cartoon': {'color':color}})
                                xyzview1.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "lightgrey"},{"hetflag": False})

                            else:
                                xyzview1.setStyle({'chain':chain},{'stick': {'color':color, "radius":  0.4}})
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
        elif st.button('Process with wiggle',key="process_wiggle"):
            st.write('')
            s=time.time()
            with st.spinner('Processing...'):
                g = algo.attachwithwiggle(protein,glycans,glycosylation_locations)
            st.balloons()
            st.success("exec time :"+str(time.time()-s)+"Seconds")
            g1 = pdb.exportPDB('output/out.pdb',pdb.to_normal(g))
            xyzview1 = py3Dmol.view()
            xyzview1.addModelsAsFrames(g1)
            for n,chain,color in zip(range(len(glycosylation_locations)+1),list("ABCDEFGH"),
                                        ["grey","#FF7B89","#8A5082","#6F5F90","#758EB7","#A5CAD2","blue","orange"]):
                            if chain=="A":
                                xyzview1.setStyle({'chain':chain},{'cartoon': {'color':color}})
                                xyzview1.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "lightgrey"},{"hetflag": False})

                            else:
                                xyzview1.setStyle({'chain':chain},{'stick': {'color':color, "radius":  0.4}})
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
    st.write(bytes_data)    
else:
    st.header('Just Another Glycosylator for Nerds')
    st.subheader("Yeahh! It's uses MD simulation results, so, obviously its better ... Duhh!")
    st.markdown("")
    st.markdown("")

