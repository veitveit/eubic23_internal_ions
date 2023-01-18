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
import plotly.figure_factory as ff
import pandas as pd
import json 




def example_plot():
    x1 = np.random.randn(200) - 2
    x2 = np.random.randn(200)
    x3 = np.random.randn(200) + 2
    hist_data = [x1, x2, x3]
    group_labels = ['Group 1', 'Group 2', 'Group 3']
    fig = ff.create_distplot(hist_data, group_labels, bin_size = [.1, .25, .5])
    return fig






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
        
    with data_import_tab:

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
    
        Spectrum_select = st.number_input("Input desired spetrum number" , 
                                          step = 1, min_value = 1, max_value = len(pandas_df))

    

    with stats_tab:   
    
        if st.button("Make some plots!", help = "desc"):
            plot = st.plotly_chart(example_plot(), use_container_width = True) 

        
        
    
        

# side bar and main page loader
def main():

    about_str = \
    """
    description
    """

    st.set_page_config(page_title = "Internal ions",
                       page_icon = ":test_tube:",
                       layout = "centered",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/.../discussions",
                                     "Report a bug": "https://github.com/.../issues",
                                     "About": about_str}
                       )

    title = st.sidebar.title("Internal ions")

    logo = st.sidebar.image("img/logo.png", caption = "Logo")

    doc = st.sidebar.markdown("desc")
    
    #pages = ("Main page", "Data import", "Spectra analysis", "Statistical analysis")
    
    #page = st.sidebar.selectbox(label = "Start the internal fragment analysis",
    #                            options = pages,
    #                            index = 0,
     #                           help = "Select a workflow that you want to run.")
                                
    
    
    contact_str = "**Contact:** [Micha Birklbauer](mailto:micha.birklbauer@gmail.com)"
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/michabirklbauer/piaweb/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)

    
    main_page()
    

        
    
    
    
def run_fragannot(spectrum_file, identifications_file): 
    """
    Funtion that runs fragannot and outputs a json file.  
    
    
    
    """
    
    json_file = "file"

    return json_file



def json_pandas(json_file): 
    """ 
    Function to read the json fragannot output file into a pandas df. 
    
    """
    
    # Opening JSON file
    f = open(json_file)

    # returns JSON object as 
    # a dictionary
    data = json.load(f)

    # Iterating through the json
    # list

    f.close()
    
    return pd.DataFrame(data)


def filter_json_file(pd_df, start_seq_length, end_seq_length, start_frag_len, end_frag_len, start_mz, end_mz, start_int, end_int): 
    """ 
    
    """
    filter_df = "file"
    
    return filter_df 



    
    
    

if __name__ == "__main__":
    main()
