#webframework
import streamlit as st

#to display dataframe 
import pandas as pd

#display image
from PIL import Image

#allows descriptor calculation
import subprocess
import os
import base64

#to calculate lipinski descriptors
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import Crippen
import numpy as np

#to load pickle file
import pickle

#to convert IUPAC names to smiles
from urllib.request import urlopen
from urllib.parse import quote
from pathlib import Path
import webbrowser

#to plot graph
import plotly.express as px
from plotly.subplots import make_subplots

# Molecular descriptor calculation
def descriptor_calculation():
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_final_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')
 
#linpinksi calculation 
def lipinski_iupac(smiles, verbose=False):
    mol_data= [] 
        
    try:
        mol=Chem.MolFromSmiles(smiles) 
        mol_data.append(mol)
    except:
        print('Invalid SMILES:', smiles)
        
       
    base_data= np.arange(1,1)
                
    desc_MolWt = Descriptors.MolWt(mol)
    desc_MolLogP = Descriptors.MolLogP(mol)
    desc_NumHDonors = Lipinski.NumHDonors(mol)
    desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
    row = np.array([desc_MolWt,
                    desc_MolLogP,
                    desc_NumHDonors,
                    desc_NumHAcceptors])   
    

    base_data=row
    
    column_names=["MW","LogP","HBD","HBA"]
    descriptors = pd.DataFrame([base_data],columns=column_names)
    
    return descriptors

def lipinski(smiles, verbose=False):
    mol_data= []    
    for element in smiles:
     
      try:
          mol=Chem.MolFromSmiles(element) 
          mol_data.append(mol)
      except:
          print('Invalid SMILES:', element)
        
       
    base_data= np.arange(1,1)
    i=0  
    for mol in mol_data:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            base_data=row
        else:
            base_data=np.vstack([base_data, row])
        i=i+1      
    
    column_names=["MW","LogP","HBD","HBA"]   
    descriptors = pd.DataFrame(data=base_data,columns=column_names)
    
    return descriptors

# Model building
def model_predict(input_data):
    load_model = pickle.load(open('bioactivity_prediction_model.pkl', 'rb'))
    prediction = load_model.predict(input_data)
    
    st.header('**Prediction output**')
    
    prediction_output = pd.Series(prediction, name='pIC50')
    
    if identifier!="":
        molecule_name = pd.Series(iupac_load_data, name='molecule_name')
        df = pd.concat([molecule_name, prediction_output], axis=1)
    elif identifier2!="":
        molecule_name = pd.Series(draw_data, name='molecule_name')
        df = pd.concat([molecule_name, prediction_output], axis=1)
    else:
        molecule_name = pd.Series(load_data[1], name='molecule_name')
        oral=[] 
        for i in load_data[0]:
            result=lipinski_pass(i)            
            if result==True:
                oral.append("Orally Bioavailable")
            else:
                oral.append("Not Orally Bioavailable")
        sm_col = pd.Series(oral, name='Oral Bioavailability')
        df = pd.concat([molecule_name, prediction_output,sm_col], axis=1)
    
    
    st.write(df)
    st.markdown(download_csv(df), unsafe_allow_html=True)
    
    df = pd.DataFrame(df)
    return df
    
# download prediction
def download_csv(df):
    csv = df.to_csv(index=False)    # if no filename is given then a string is returned
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download</a>'
    return href 

#convert iupac name to smiles
def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'

#Returns the octanol-water partition coefficient given a molecule SMILES string
def log_partition_coefficient(smiles):
       
    mol = Chem.MolFromSmiles(smiles)        
    return Crippen.MolLogP(mol)
    
def lipinski_trial(smiles,verbose=False):
    '''    
    Lipinski's rule:
    Hydrogen bond donors <= 5
    Hydrogen bond acceptors <= 10
    Molecular weight < 500 daltons
    logP < 5
    '''
    lip_pass = []
    lip_fail = []
    
    mol = Chem.MolFromSmiles(smiles)
    
    num_hdonors = Lipinski.NumHDonors(mol)
    num_hacceptors = Lipinski.NumHAcceptors(mol)
    mol_weight = Descriptors.MolWt(mol)
    mol_logp = Crippen.MolLogP(mol)
    
    lip_fail = []
    
    if num_hdonors > 5:
        lip_fail.append('Over 5 H-bond donors')
    else:
        lip_pass.append('Found %s H-bond donors')
        
    if num_hacceptors > 10:
        lip_fail.append('Over 10 H-bond acceptors')
    else:
        lip_pass.append('Found %s H-bond acceptors')
        
    if mol_weight >= 500:
        lip_fail.append('Molecular weight over 500')
    else:
        lip_pass.append('Molecular weight: %s')
        
    if mol_logp >= 5:
        lip_fail.append('Log partition coefficient over 5')
    else:
        lip_pass.append('Log partition coefficient: %s')
    
    return lip_pass, lip_fail
    
def lipinski_pass(smiles):
    lip_pass, lip_fail = lipinski_trial(smiles)
    if lip_fail:
        return False
    else:
        return True    

image = Image.open('drug.jpg')
st.image(image, use_column_width=True)

st.markdown("""
# Bioactivity Prediction Application: Acetylcholinesterase

***Predict the bioactivity towards inhibiting Acetylcholinesterase responsible for diseases like Parkinson's, schizophrenia and Alzheimer's.***

""")


with st.sidebar.header('Enter input data'):
    uploaded_file = st.sidebar.file_uploader("1. Upload your input file with smiles notation of molecules", type=['txt'])
    

with st.sidebar.header('2. Enter IUPAC name of molecule'):
    iupac_name = st.text_input("2. Enter IUPAC name of molecule")
    identifier  = iupac_name
    
with st.sidebar.header('3. Draw molecule by using below link'):
    st.write('3.Draw molecule by using below link')
    
with st.sidebar.header('3. Draw molecule by using below link'):
    url = 'http://localhost:5006/jsme'

    if st.button('Draw'):
         webbrowser.open_new_tab(url)
       
with st.sidebar.header('Paste the smile notation'):
    draw_data = st.text_input("Paste the smile notation")
    identifier2  = draw_data

flag=1
if st.sidebar.button('Predict'):

    if identifier!="":
        iupac_load_data = CIRconvert(identifier)
        if iupac_load_data=='Did not work':
            st.warning('Did not work. Please enter another name.')
            flag=0
        else:
            path_to_file = 'molecule.smi'
            path = Path(path_to_file)

            if path.is_file():    
                os.remove('molecule.smi')
            else:
                with open("molecule.smi","x") as file:
                    file.write(iupac_load_data + "\n")
            st.header('**Input data**')
            st.write(identifier, '  =>  ',iupac_load_data)
    
    elif identifier2!="":
        path_to_file = 'molecule.smi'
        path = Path(path_to_file)

        if path.is_file():    
            os.remove('molecule.smi')
        else:
            with open("molecule.smi","x") as file:
                file.write(draw_data + "\n")
        st.header('**Input data**')
        st.write(draw_data)
    
    else:    
        load_data = pd.read_table(uploaded_file, sep=' ', header=None)
        temp = load_data
        load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

        st.header('**Input data**')
        st.write(load_data) 
 
    #calculating molecular descriptors
    if flag==1:
        with st.spinner("Calculating..."):
            descriptor_calculation()

        # Read in calculated descriptors and display the dataframe
        st.header('**Calculated molecular descriptors**')
        desc = pd.read_csv('descriptors_final_output.csv')
        
        if identifier!="":
            temp = lipinski_iupac(iupac_load_data)
        elif identifier2!="":
            temp = lipinski_iupac(draw_data)
        else:
            temp = lipinski(temp[0])
        
        desc_final = pd.concat([desc,temp], axis=1)
        
        st.write(desc_final)
        st.write(desc_final.shape)

        st.header('**Selected molecular descriptors**')
        Xlist = list(pd.read_csv('descriptor_final_list.csv').columns)
        desc_subset = desc_final[Xlist]
        
        st.write(desc_subset)
        st.write(desc_subset.shape)

        # Apply trained model to make prediction 
        with st.spinner("Model Predicting..."):        
            graph = model_predict(desc_subset)
        
        #Check for oral bioavailabiliy
        if identifier!="" or identifier2!="":
            if identifier!="":
                result=lipinski_pass(iupac_load_data)
            else:
                result=lipinski_pass(draw_data)
            if result==True:
                new_title = '<p style="color:Green; font-size: 42px;">*Orally bioavailable*</p>'
                st.markdown(new_title, unsafe_allow_html=True)
            else:
                new_title = '<p style="color:Red; font-size: 42px;">*Not orally bioavailable*</p>'
                st.markdown(new_title, unsafe_allow_html=True)            

        #plotting graph and finding most compatible drug
        if  identifier=="" and identifier2=="":
            st.header('**pIC50 Graph**')        
            with st.spinner("Plotting graph..."): 
                graph_output = px.bar(
                    graph, 
                    x='molecule_name',
                    y='pIC50',
                    color = "molecule_name")
                st.plotly_chart(graph_output)
            st.header('**Most compatible molecule to be used as drug**')
            column_p=graph["pIC50"]
            column_o=graph["Oral Bioavailability"]
            i=0
            index=0
            max_val=0
            while i<len(graph):
                if column_p[i]>max_val and column_o[i]=="Orally Bioavailable":
                    max_val = column_p[i]
                    index=i
                i=i+1
            comp_mol=graph["molecule_name"]            
            st.write(comp_mol[index],'  => ',max_val)        
    
else:
    st.info('Upload input data in the sidebar for prediction')
       
#2,5,5-trimethyl-2-hexene