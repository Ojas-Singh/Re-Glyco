import streamlit as st
import pandas as pd
import numpy as np
from lib import prep,algo,pdb
import py3Dmol
from stmol import showmol


st.set_page_config(layout="wide")
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
        style = st.sidebar.selectbox('style',['line','cross','stick','sphere','cartoon','clicksphere'])
        xyzview = py3Dmol.view()
        xyzview.addModelsAsFrames(system)
        xyzview.setStyle({style:{'color':'spectrum'}})
        xyzview.setBackgroundColor('#0e1117')
        xyzview.zoomTo()
        showmol(xyzview,height=400,width=800)
        expander = st.expander("See glycosylation locations!")
        expander.write(out["glycosylations"])
        glycosylation_locations = out["glycosylations"]
        options = st.multiselect(
                    'What Glycans to attach?',
                    ['bisecting','Man3','Man5','Man9'],
                    ['bisecting'])
        if st.button('Process'):
            st.write('')
            g = algo.attach(protein,options,glycosylation_locations)
            g1 = pdb.exportPDB('output/out.pdb',pdb.to_normal(g))
            print(g)
            print("ok")
            xyzview1 = py3Dmol.view()
            xyzview1.addModelsAsFrames(g1)
            xyzview1.setStyle({'stick':{'color':'spectrum'}})
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
    print("Enter UniProtID to procede!")

