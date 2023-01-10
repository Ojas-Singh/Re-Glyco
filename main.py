import streamlit as st
import pandas as pd
import numpy as np
from lib import prep,algo,pdb
import py3Dmol
from stmol import showmol
from PIL import Image


st.set_page_config(layout="wide")
image = Image.open('logo.png')
st.sidebar.image(image, caption='')
st.sidebar.title('Re-Glyco')
uni_id = st.sidebar.text_input(
        "Enter UNIPROT ID of protien...","O15552"
    )

@st.cache
def fetch(uni_id):
    fold = prep.download_and_prepare_alphafoldDB_model(uni_id,"output/temp/")
    out= prep.query_uniprot_for_glycosylation_locations(uni_id)
    return fold,out


if uni_id != "":
    fold, out = fetch(uni_id)
    with open(fold) as ifile:
        system = "".join([x for x in ifile])
        btn = st.download_button(
                label="Download AlphaFold Structure",
                data=system,
                file_name=uni_id+".pdb",
                mime='text/csv'
            )
        protein = pdb.parse(fold)
        style = st.sidebar.selectbox('style',['cartoon','stick'])
        xyzview = py3Dmol.view()
        xyzview.addModelsAsFrames(system)
        xyzview.setStyle({style:{'color':'spectrum'}})
        xyzview.setBackgroundColor('#0e1117')
        xyzview.zoomTo()
        showmol(xyzview,height=400,width=800)
        expander = st.expander("See glycosylation locations!")
        expander.write(out["glycosylations"])
        glycosylation_locations = out["glycosylations"]
        options = st.selectbox(
                    'What Glycans to attach?',
<<<<<<< HEAD
                    ('bisecting', 'a2',"a2g2","a3g3","m5","m6_1","m6_2","m6_3","m7_1","m7_2","m7_3","m7_4","m8_1","m8_2","m8_3","m9"))
        attempt = st.number_input('Number of iterations :',min_value=100, max_value=5000, value=1000, step=500, format=None)
        if st.button('Process'):
=======
                    ('bisecting', 'man'))
        # attempt = st.number_input('Number of iterations :',min_value=100, max_value=50000, value=1000, step=500, format=None)
        if st.button('Process without wiggle'):
>>>>>>> 0da8f9d1bd5599e1c63f9b0b51c4ae01c8f606a8
            st.write('')
            g = algo.attach(protein,options,glycosylation_locations)
            g1 = pdb.exportPDB('output/out.pdb',pdb.to_normal(g))
            # print(g)
            print("ok")
            xyzview1 = py3Dmol.view()
            xyzview1.addModelsAsFrames(g1)
            # import random
            # colors = {i: "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(120)}
            # xyzview1.setStyle({'stick':{'colorscheme':{'prop':'chain','map':colors}}})
            for n,chain,color in zip(range(len(glycosylation_locations)+1),list("ABCDEFGH"),
                                        ["lime","cyan","magenta","yellow","salmon","white","blue","orange"]):
                            if chain=="A":
                                xyzview1.setStyle({'chain':chain},{'cartoon': {'color':color}})
                            else:
                                xyzview1.setStyle({'chain':chain},{'stick': {'color':color}})


            # xyzview1.setStyle({'stick':{'color':'spectrum'}})
            xyzview1.setBackgroundColor('#0e1117')
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
        elif st.button('Process with wiggle'):
            st.write('')
            g = algo.attachwithwiggle(protein,options,glycosylation_locations)
            g1 = pdb.exportPDB('output/out.pdb',pdb.to_normal(g))
            # print(g)
            print("ok")
            xyzview1 = py3Dmol.view()
            xyzview1.addModelsAsFrames(g1)
            # import random
            # colors = {i: "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(120)}
            # xyzview1.setStyle({'stick':{'colorscheme':{'prop':'chain','map':colors}}})
            for n,chain,color in zip(range(len(glycosylation_locations)+1),list("ABCDEFGH"),
                                        ["lime","cyan","magenta","yellow","salmon","white","blue","orange"]):
                            if chain=="A":
                                xyzview1.setStyle({'chain':chain},{'cartoon': {'color':color}})
                            else:
                                xyzview1.setStyle({'chain':chain},{'stick': {'color':color}})


            # xyzview1.setStyle({'stick':{'color':'spectrum'}})
            xyzview1.setBackgroundColor('#0e1117')
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
        
        

else:
    print("Enter UniProtID to proceed!")

