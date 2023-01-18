import streamlit as st
import pandas as pd
import numpy as np
from lib import prep,algo,pdb
import py3Dmol
from stmol import showmol
from PIL import Image
import matplotlib.pyplot as plt
import plotly.express as px
from sklearn.cluster import KMeans
import scipy.stats as stg
import plotly.graph_objects as go

import time
# st.set_page_config(layout="wide")
# image = Image.open('logo.png')
# st.sidebar.image(image, caption='')
st.title('Glycan Conformation Sampler for Nerds')
st.header('')
st.subheader("")
name="bisecting"
fold="data/bisecting.pdb"
df = pd.read_csv('output/'+name+'torsions.csv')
pca_df = pd.read_csv('output/'+name+'pca.csv')
with open(fold) as ifile:
    system = "".join([x for x in ifile])
    tab1, tab2, tab3 = st.tabs(["Structure", "Clusters", "Sampler"])
    with tab1:
        col1, col2 = st.columns(2)
        with col1:
            pass
        protein = pdb.parse(fold)
        xyzview = py3Dmol.view()
        xyzview.addModelsAsFrames(system)
        xyzview.setStyle({'stick':{'color':'spectrum'}})
        xyzview.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "lightgrey"},{"hetflag": False})

        xyzview.setBackgroundColor('#FFFFFF')
        xyzview.zoomTo()
        showmol(xyzview,height=800,width=900)

    with tab2:
        fig0 = px.scatter(
            pca_df,
            x="X",
            y="Y",
            color="i",
            # color_continuous_scale="reds",
        )
        fig0.update_traces(marker=dict(size=2,),
                  selector=dict(mode='markers'))
        st.plotly_chart(fig0, theme="streamlit", use_conatiner_width=True)
        n_clusters = st.slider('How many clusters?', 0, 50, 10)
        clustering = KMeans(n_clusters).fit(pca_df[['X','Y']])
        clustering_labels= clustering.labels_
        pca_df.insert(1,"cluster",clustering_labels,False)
        df.insert(1,"cluster",clustering_labels,False)
        
        fig1 = px.scatter(
            pca_df,
            x="X",
            y="Y",
            color="cluster",
            # color_continuous_scale="reds",
        )
        fig1.update_traces(marker=dict(size=2,),
                  selector=dict(mode='markers'))
        st.plotly_chart(fig1, theme="streamlit", use_conatiner_width=True)
        popp=[]
        for i in range(n_clusters):
            df0 = pca_df.loc[df["cluster"] ==int(i)]
            o=[]
            pp=clustering.cluster_centers_[i]
            for j in range(len(df0.iloc[:,0])):
                o.append([np.linalg.norm(np.asarray(pp)-[df0["X"].iloc[j],df0["Y"].iloc[j]]),df0["i"].iloc[j]])
            o.sort()
            popp.append(o[0][1])
        st.write(popp)

        xax = st.selectbox(
    'How would you like to be contacted?',
    (list(df.columns.values)),key="x")
        yax = st.selectbox(
        'How would you like to be contacted?',
        (list(df.columns.values)),key="y")
        t = np.linspace(-1, 1.2, 2000)
        psi = df[xax]
        phi = df[yax]
        x=psi
        y=phi
        fig = go.Figure(go.Histogram2dContour(
                x = x,
                y = y,
                colorscale = 'Blues'
        ))

        st.plotly_chart(fig, theme="streamlit", use_container_width=True)
        fig2 = px.scatter(
            df,
            x=xax,
            y=yax,
            color="cluster",
            # color_continuous_scale="reds",
        )
        fig2.update_traces(marker=dict(size=2,),
                  selector=dict(mode='markers'))
        st.plotly_chart(fig2, theme="streamlit", use_conatiner_width=True)
        

    with tab3:
        # clustering = KMeans(50).fit(pca_df[['X','Y']])
        # clustering_labels= clustering.labels_
        # pca_df.insert(1,"cluster",clustering_labels,False)
        # df.insert(1,"cluster",clustering_labels,False)
        
        fig1 = px.scatter(
            pca_df,
            x="X",
            y="Y",
            color="cluster",
            # color_continuous_scale="reds",
        )
        fig1.update_traces(marker=dict(size=2,),
                  selector=dict(mode='markers'))
        st.plotly_chart(fig1, theme="streamlit", use_conatiner_width=True)
        cluster1= st.selectbox(
    'From cluster?',
    (list(range(50))),key="clu")
        torsionrange = st.slider('Random torsion range?', 0, 10, 3)
        import os
        try:
            os.remove('output/wig.pdb') 
        except:
            pass
        if st.button('Process',key="process"):
            G = pdb.parse("output/cluster/"+str(cluster1)+".pdb")
            G= pdb.to_DF(G)
            loaded = np.load('data/bisecting.npz',allow_pickle=True)
            Garr = G[['X','Y','Z']].to_numpy(dtype=float)
            tormeta = loaded["b"]

            # torsions = loaded["c"]

            torsionpoints = loaded["d"]
            torsionparts  = loaded["f"]
            torsionparts = np.asarray(torsionparts)
            torsionpoints= np.asarray(torsionpoints)
            molecules = []
            for idx in range(20):
                Garr1 = algo.Garrfromtorsiondemo(Garr,torsionpoints,torsionrange,torsionparts)
                Gn =  pd.DataFrame(Garr1, columns = ['X','Y','Z'])
                G.update(Gn)
                g1 = pdb.exportPDBmulti('output/wig.pdb',pdb.to_normal(G),idx)
            with open('output/wig.pdb') as ifile:
                systemx = "".join([x for x in ifile])
                xyzview1 = py3Dmol.view()
                xyzview1.addModelsAsFrames(systemx)
                xyzview1.setStyle({'stick':{'color':'spectrum'}})
                # xyzview1.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "lightgrey"},{"hetflag": False})

                xyzview1.setBackgroundColor('#FFFFFF') 
                xyzview1.zoomTo()
                xyzview1.animate({'loop': "forward"})
                xyzview1.show()
                showmol(xyzview1,height=800,width=900)