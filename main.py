import streamlit as st
from lib import prep,algo,pdb,molview
import time,os,config
from PIL import Image
from io import StringIO
from annotated_text import annotated_text
import streamlit.components.v1 as components
import base64


def get_glycan_list(suffix):
    directory = config.data_dir
    folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder)) and folder.endswith(suffix)]
    return folders


st.set_page_config(page_title="Re-GLyco", page_icon='resource/icon.png', layout="wide", initial_sidebar_state="auto", menu_items=None)

hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 
image = Image.open('static/resource/logo.png')
st.sidebar.image(image, caption='')
uniprot=""
uni_id = st.sidebar.text_input(
        "Enter UNIPROT ID of protein...",placeholder ="P29016 etc"
    )

st.sidebar.write("Or upload PDB manually.")
uploaded_file = st.sidebar.file_uploader("")


if not uni_id=="" and uploaded_file is None:
    fold, out = prep.fetch(uni_id)
    with open(fold) as ifile:
        system = "".join([x for x in ifile])
        confidence =pdb.get_confidence(system)
        tab1, tab2 = st.tabs(["Structure", "glycosylations"])
        glycosylation_locations = out["glycosylations"]
        with tab1:
            col1, col2 = st.columns(2)
            with col1:
                molview.show3d(fold,system,confidence,glycosylation_locations)
                protein = pdb.parse(fold)
                
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
            st.write(glycosylation_locations)
        glycans=[1 for x in range(len(glycosylation_locations))]
        glycosylation_locations_N=[]
        for i in range(len(glycosylation_locations)):
            if not (glycosylation_locations[i]["description"].startswith('N-linked') or glycosylation_locations[i]["description"].startswith('O-linked') or glycosylation_locations[i]["description"].startswith('C-linked')):
                st.warning(f"Only N & O glycosylation supported, location : {glycosylation_locations[i]['begin']}  is {glycosylation_locations[i]['description']}")
            else:
                glycosylation_locations_N.append(glycosylation_locations[i])
        dirlist_N = ["None",] + get_glycan_list("Man(b1-4)GlcNAc(b1-4)GlcNAc")
        dirlist_O = ["None",] + get_glycan_list("GalNAc") 
        dirlist_C = ["None",] + get_glycan_list("Man")
        for i in range(len(glycosylation_locations_N)): 
            df = pdb.to_DF(protein)
            resname = df.loc[df['ResId'] == int(glycosylation_locations_N[i]["begin"]), 'ResName'].iloc[0]
            if resname in config.N_linked["Res"]:
                options = st.selectbox(
                            f'Which glycans to attach on location : {resname}{glycosylation_locations_N[i]["begin"]} ?',
                            (dirlist_N),key=str(i))
                glycans[i]=options
            elif resname in config.O_linked["Res"]: 
                options = st.selectbox(
                            f'Which glycans to attach on location : {resname}{glycosylation_locations_N[i]["begin"]} ?',   
                            (dirlist_O),key=str(i))
                glycans[i]=options
            elif resname in config.C_linked["Res"]:
                options = st.selectbox(
                            f'Which glycans to attach on location : {resname}{glycosylation_locations_N[i]["begin"]} ?',   
                            (dirlist_C),key=str(i))
                glycans[i]=options
        if st.button('Process',key="process") and not all(element is "None" for element in glycans):
            s=time.time()
            with st.spinner('Processing...'):
                g,clash = algo.attach(protein,glycans,glycosylation_locations_N)
            if clash:
                clashh=True
                st.warning('Clash detected, rerun or the spot is not glycosylable, [low confidence region near spot.]  ')
            st.balloons()
            st.success("exec time : "+ str(int(time.time()-s)) +" seconds")
            g1 = pdb.exportPDB('/mnt/database/server_dir/temp/out.pdb',pdb.to_normal(g))
            # components.iframe("https://healoor.me/litemol/index.html?pdbUrl=https://glycoshape.healoor.me/temp/out.pdb",height=600)
            molview.show3doutput(g1,glycosylation_locations,confidence)
            with open('/mnt/database/server_dir/temp/out.pdb') as ofile:
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
    confidence= pdb.get_confidence(string_data)
    protein = pdb.parse("output/temp/custom.pdb")
    molview.show3dbasic(string_data)
    protein_df= pdb.to_DF(protein)
    seqeunce,shift = pdb.get_sequence_from_pdb("output/temp/custom.pdb")
    st.code(seqeunce)
    pridicted_spots = pdb.find_glycosylation_spots(seqeunce,shift)
    # spots= protein_df.loc[(protein_df['ResName']=="ASN") |(protein_df['ResName']=="THR") | (protein_df['ResName']=="TYR")|(protein_df['ResName']=="TRP") | (protein_df['ResName']=="SER") ,['ResId']].iloc[:]['ResId'].tolist()
    spots= protein_df.loc[(protein_df['ResName']=="THR") | (protein_df['ResName']=="TYR")|(protein_df['ResName']=="TRP") | (protein_df['ResName']=="SER") ,['ResId']].iloc[:]['ResId'].tolist()

    def unique(spots):
        unique_list = []
        for x in spots:
            if x not in unique_list:
                unique_list.append(x)
        return unique_list
    spots = pridicted_spots + unique(spots) 
    glycosylation_locations = st.multiselect(
    'Select Glycosylation Locations',
    spots,
    )
    glycans=[1 for x in range(len(glycosylation_locations))]
    df = pdb.to_DF(protein)
    dirlist_N = ["None",] + get_glycan_list("Man(b1-4)GlcNAc(b1-4)GlcNAc")
    dirlist_O = ["None",] + get_glycan_list("GalNAc") 
    dirlist_C = ["None",] + get_glycan_list("Man")
    for i in range(len(glycosylation_locations)): 
        df = pdb.to_DF(protein)
        resname = df.loc[df['ResId'] == int(glycosylation_locations[i]), 'ResName'].iloc[0]
        if resname in config.N_linked["Res"]:
            options = st.selectbox(
                        f'Which glycans to attach on location : {resname}{glycosylation_locations[i]}  ?',
                        (dirlist_N),key=str(i))
            glycans[i]=options
        elif resname in config.O_linked["Res"]: 
            options = st.selectbox(
                        f'Which glycans to attach on location : {resname}{glycosylation_locations[i]} ?',   
                        (dirlist_O),key=str(i))
            glycans[i]=options
        elif resname in config.C_linked["Res"]:
            options = st.selectbox(
                        f'Which glycans to attach on location : {resname}{glycosylation_locations[i]} ?',   
                        (dirlist_C),key=str(i))
            glycans[i]=options
    if st.button('Process',key="process") and not all(element is "None" for element in glycans):
            s=time.time()
            with st.spinner('Processing...'):
                g,clash = algo.attach(protein,glycans,glycosylation_locations)
            if clash:
                clashh=True
                st.warning('Clash detected, rerun or the spot is not glycosylable, [low confidence region near spot.]  ')
            st.balloons()
            st.success("exec time :"+str(time.time()-s)+"Seconds")
            g1 = pdb.exportPDB('/mnt/database/server_dir/temp/out.pdb',pdb.to_normal(g))
            # components.iframe("https://healoor.me/litemol/index.html?pdbUrl=https://glycoshape.healoor.me/temp/out.pdb",height=600)
            molview.show3doutput(g1,glycosylation_locations,confidence)
            with open('/mnt/database/server_dir/temp/out.pdb') as ofile:
                system = "".join([x for x in ofile])
                btn = st.download_button(
                    label="Download Re-Glyco Structure",
                    data=system,
                    file_name=uni_id+"_glycosylated.pdb",
                    mime='text/csv'
                )
else:
    st.title("Welcome to re-glyco!")
    # st.header('Under Construction!')
    st.markdown("""
    Re-glyco is a powerful tool designed to restore the missing glycosylation in AlphaFold structures or user-uploaded protein structures.
    
    To get started, upload your protein structure file or choose a pre-existing AlphaFold structure, and let re-glyco do the rest!

    here are some example UniProt IDs to get you started:

    O15552, P29016, Q9BXJ4, P27918, B0YJ81
    """)
    
