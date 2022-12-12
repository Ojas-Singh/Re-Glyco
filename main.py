import streamlit as st
import pandas 
import numpy
from lib import prep,algo,pdb
import py3Dmol
from stmol import showmol


st.set_page_config(layout="wide")
st.sidebar.title('Re-Glyco')
uni_id = st.sidebar.text_input(
        "Enter UNIPROT ID of protien...","O15552"
    )

if uni_id != "":
    fold = prep.download_and_prepare_alphafoldDB_model(uni_id,"output/temp/")
    out= prep.query_uniprot_for_glycosylation_locations(uni_id)
        
    with open(fold) as ifile:
        system = "".join([x for x in ifile])
        btn = st.download_button(
                label="Download AlphaFold Structure",
                data=system,
                file_name=uni_id+".pdb",
                mime='text/csv'
            )
        protien = pdb.fullparse(fold)
        style = st.sidebar.selectbox('style',['line','cross','stick','sphere','cartoon','clicksphere'])
        xyzview = py3Dmol.view()
        xyzview.addModelsAsFrames(system)
        xyzview.setStyle({style:{'color':'spectrum'}})
        xyzview.setBackgroundColor('#0e1117')
        xyzview.zoomTo()
        showmol(xyzview,height=400,width=800)
        expander = st.expander("See glycosylation locations!")
        expander.write(out)
        st.write()
else:
    print("Enter UniProtID to procede!")

