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
st.sidebar.title('Wiggle')

