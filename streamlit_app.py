#!/usr/bin/env python3

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

import numpy as np
import streamlit as st

import pandas as pd
import json 

from json_to_dataframes import json_to_dataframes


from stats_page import * 
from spectra_page import * 
from filtering_files import * 
    




# main page content
def main_page(): 
    
    main_tab, data_import_tab, spectra_tab, stats_tab = st.tabs(["Main page", "Data import", 
                                                                  "Spectra analysis", "Statistical analysis"])

    with main_tab:
        title = st.title("Internal ions")

        text_1 = st.markdown("General description")   

        # Decripbe project, goal and members 

        text_1 = st.markdown("Internal ion description ") 
        
        # Add forumula for nomenclatur and mass 
        formula = st.latex(r'''
                    a + ar + a r^2 + a r^3 + \cdots + a r^{n-1} =
                    \sum_{k=0}^{n-1} ar^k =
                    a \left(\frac{1-r^{n}}{1-r}\right)
                    ''')

        
        
    with data_import_tab:
        
        title = st.title("Data Import")

        spectrum_file = st.file_uploader("Upload spectrum file:",
                                         type = ["raw"],
                                         help = "desc"
                                        )
    
        identifications_file = st.file_uploader("Upload identification file:",
                                                type = ["mzId"],
                                                help = "desc"
                                               )

   
        # Running fragannot 
        #json_data = run_fragannot(spectrum_file, identifications_file)  

        # Change this later 
        test_FILE = "/Users/hmt128/Work/eubic/data/data.json" 

        # converting json into pandas 
        pandas_df = json_pandas(test_FILE)  # change to json_data 

        data_text = st.markdown("Your data has arrived")   

        st.table(pandas_df.head(10))

        data_text = st.markdown("Select filtering options below")  

        start_seq_length, end_seq_length = st.select_slider("Peptide sequence lenght", 
                                            options= range(0,501), 
                                            value = (10,50))

        start_frag_len, end_frag_len = st.select_slider("Fragment ion sequence lenght", 
                                            options= range(0,101), 
                                            value = (10,50))

        start_mz, end_mz = st.select_slider("Select m/z range", 
                                            options= range(0,1001), 
                                            value = (20,50))

        start_int, end_int = st.select_slider("Select intensity range", 
                                              options= range(0,100001), 
                                              value =(20,10000))

        #seq_length = st.number_input("Identification score range")
        data_text = st.markdown("Select ions to look for below") 
        A_ion = st.checkbox("A ions")
        B_ion = st.checkbox("B ions")
        C_ion = st.checkbox("C ions")
        X_ion = st.checkbox("X ions")
        Y_ion = st.checkbox("Y ions")
        Z_ion = st.checkbox("Z ions")
        
        # running filtering function 
    
        filt_df = filter_json_file(pandas_df, start_seq_length, 
                               end_seq_length, start_frag_len, end_frag_len, 
                                   start_mz, end_mz, start_int, end_int) 
    
    
    with spectra_tab:
        
        title = st.title("Spectra view")
        
        data_text = st.markdown("Desc")
        
    
        Spectrum_select = st.number_input("Input desired spetrum number" , 
                                          step = 1, min_value = 1, max_value = len(pandas_df)) 
        
        # Make a function thart shows one spectra at a time 
        # select one spectra from dataframe to pass into this fucntion
        
        show_specta(Spectrum_select)

    

    with stats_tab:
        
        title = st.title("Statistics view")
        
        data_text = st.markdown("Desc")
    
        title = st.title("First plot")
        plot = st.plotly_chart(common_type_hist(test_FILE), use_container_width = True) 
            
        title = st.title("Second plot")
        plot = st.plotly_chart(example_plot(), use_container_width = True) 

        title = st.title("Second plot")
        plot = st.plotly_chart(example_plot(), use_container_width = True) 
        
    


# side bar and main page loader
def main():

    about_str = \
    """
    description
    """

    st.set_page_config(page_title = "Internal ions",
                       page_icon = "random", # ":test_tube:"
                       layout = "wide",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/.../discussions",
                                     "Report a bug": "https://github.com/.../issues",
                                     "About": about_str},  # 

                       )

    title = st.sidebar.title("EuBic 2023 Hackathon")

    logo = st.sidebar.image("img/logo.png", caption = "Logo")

    doc = st.sidebar.markdown("desc")
   
                               
    
    contact_str = "**Contact:** [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk), [Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk), [Kristian Boje Nielsen](krisn16@student.sdu.dk), [Louise Buur](louise.buur@fh-hagenberg.at), [Micha Birklbauer](mailto:micha.birklbauer@gmail.com), [Mohieddin Jafari](mohieddin.jafari@helsinki.fi), [Veit](veits@bmb.sdu.dk), [Vladimir Gorshkov](homer2k@gmail.com), [Zoltan Udvardy](zoltan.udvardy.ipbs@gmail.com) "
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/michabirklbauer/piaweb/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)

    
    main_page()
    

        
    
    
    
def run_fragannot(spectrum_file, identifications_file): 
    """
    Funtion that runs fragannot and outputs a json file.  
    
    
    
    """
    import subprocess 
    
    
    # something like this 

    subprocess.check_call("fragannot -i " + identifications_file + " -s " + spectrum_file, shell=True) 
    
    # will save files somewhere? 
    # Path to file 
    
    json_file = "file"

    return json_file



def json_pandas(json_file): 
    """ 
    Function to read the json fragannot output file into a pandas df. 
    
    """
    

    dataframes = json_to_dataframes(json_file)

    fragments_dataframe = dataframes[0]
    spectra_dataframe = dataframes[1]

    #fragments_dataframe.head()
    
    return fragments_dataframe





    
    
    

if __name__ == "__main__":
    main()
