#!/usr/bin/env python3

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

import numpy as np
import pandas as pd


import streamlit as st

from json_to_dataframes import json_to_dataframes
from stats_page import * 
from spectra_page import * 
from filtering_files import * 




# main page content
def main_page(): 
    
    main_tab, data_import_tab, stats_tab = st.tabs(["Main page", "Data import", 
                                                                  "Statistical analysis"]) # "Spectra analysis",  spectra_tab, 

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

        
        #spectrum_file = st.file_uploader("Upload spectrum file:",
        #                                 type = ["raw"],
        #                                 help = "desc"
        #                                )
    
        #identifications_file = st.file_uploader("Upload identification file:",
        #                                        type = ["mzId"],
        #                                        help = "desc"
        #                                       )


        # Change this later 
        #test_FILE_name = "/Users/hmt128/Work/eubic/data/data.json" 
        
        data_text = st.markdown("Upload and filter json files with internal fragment ion matches from fragannot here. https://github.com/arthur-grimaud/fragannot")  
        
        test_FILE = st.file_uploader("Upload json file:",
                                                type = ["json"],
                                                help = "Process you mzid and mgf file to fins internal fragment ions using fragannot on: https://github.com/arthur-grimaud/fragannot"
                                               )
        
        data_text = st.markdown("Select filtering options below")  

        start_seq_length, end_seq_length = st.select_slider("Peptide sequence lenght", 
                                            options= range(0,501), 
                                            value = (0,500))

        start_frag_len, end_frag_len = st.select_slider("Fragment ion sequence lenght", 
                                            options= range(0,101), 
                                            value = (0,100))

        start_mz, end_mz = st.select_slider("Select m/z range", 
                                            options= range(0,1001), 
                                            value = (0,1000))

        start_int, end_int = st.select_slider("Select intensity range", 
                                              options= range(0,100001), 
                                              value =(0,100000))

        #seq_length = st.number_input("Identification score range")
        data_text = st.markdown("Select which ions to analyse")  
        
        # List of booleans 
        
        check_col1, check_col2 = st.columns(2) 
        
        
        with check_col1: 
            n_ion = st.checkbox("n ions", value = True)
            A_ion = st.checkbox("A ions", value = True)
            B_ion = st.checkbox("B ions", value = True)
            C_ion = st.checkbox("C ions", value = True)
            X_ion = st.checkbox("X ions", value = True)
            Y_ion = st.checkbox("Y ions", value = True)
            Z_ion = st.checkbox("Z ions", value = True)
        
        with check_col2: 
            
            cdot_ion = st.checkbox("Cdot ions", value = True)
            cm1_ion = st.checkbox("C-1 ions", value = True)
            c1_ion = st.checkbox("C+1 ions", value = True)
            zdot_ion = st.checkbox("Zdot ions", value = True)
            z1_ion = st.checkbox("Z+1 ions", value = True)
            z2_ion = st.checkbox("Z+2 ions", value = True)
            Z3_ion = st.checkbox("Z+3 ions", value = True) 
            
        
        # ADD ALL IONS
        ion_filter_param = [A_ion, B_ion, C_ion, X_ion, Y_ion, Z_ion]
        
        
        if test_FILE is not None:
            
            dataframes = json_to_dataframes(test_FILE, is_file = True) 

               
            filt_dfs = filtering_files(dataframes, start_seq_length, end_seq_length, start_frag_len, end_frag_len, start_mz, end_mz, start_int, end_int, ion_filter_param) 
                


            # showing data in app
            data_text = st.markdown("Your data has arrived!")
            st.table(filt_dfs[0].head(10))
            
            
            with stats_tab:
        
                title = st.title("Statistics view")

                data_text = st.markdown("Desc")
                
                col1, col2 = st.columns(2)
                
                with col1: 
                    title = st.header("Commmon type Histogram")
                    plot = st.plotly_chart(common_type_hist(filt_dfs[0]), use_container_width = True) 
                    
                with col2: 
                    title = st.header("Common types Pie chart ")
                    plot = st.plotly_chart(common_type_pie(filt_dfs[0]), use_container_width = True) 
                    
                    
                col11, col22 = st.columns(2)
                
                with col11: 
                    title = st.header("Log Intensities Distribution")
                    plot = st.plotly_chart(log_ion_intens_dist(filt_dfs[0]), use_container_width = True)
                    
                with col22: 
                    # relative intensity to total intensity distribution of different ions 
                    title = st.header("Relative Log Intensities")
                    plot = st.plotly_chart(log_ion_intens_ridge(filt_dfs[0]), use_container_width = True)
                    
                    
                col11, col22 = st.columns(2)
                
                with col11: 
                    title = st.header("Log Intensities Distribution")
                    plot = st.plotly_chart(rel_ion_intens_prop(filt_dfs[0]), use_container_width = True)
                    
                with col22: 
                    # relative intensity to total intensity distribution of different ions 
                    title = st.header("Relative Log Intensities")
                    plot = st.plotly_chart(rel_ion_intens_prop_ridge(filt_dfs[0]), use_container_width = True)
                   
                
                title = st.title("Histogram of mz per ion type")
                plot = st.plotly_chart(mz_dist_ion_type(filt_dfs[0]), use_container_width = True) 
                

                col11, col22 = st.columns(2)
                
                with col11: 
                    title = st.header("Log Intensities Distribution")
                    plot = st.plotly_chart(rel_ion_intens_perc(filt_dfs[0]), use_container_width = True)
                    
                with col22: 
                    # relative intensity to total intensity distribution of different ions 
                    title = st.header("Relative Log Intensities")
                    plot = st.plotly_chart(rel_ion_intens_ridge(filt_dfs[0]), use_container_width = True)
                
                title = st.title("Per spectra")
           
                
                col111, col222 = st.columns(2)
                
                with col111: 
                    title = st.header("per_spec_ion_type")
                    plot = st.plotly_chart(per_spec_ion_type(filt_dfs[1]), use_container_width = True)
                    
                with col222: 
                    title = st.header("Log Intensities")
                    plot = st.plotly_chart(per_spec_ion_intens(filt_dfs[1]), use_container_width = True)
                
            

        #elif test_FILE is None: 
            #data_text = st.markdown("Please upload a file for analysis!")
            #None 
        

    
    #with spectra_tab:
        
        # WIP
        
        #title = st.title("Spectra view")
        
        #data_text = st.markdown("Desc")
        
    
        #Spectrum_select = st.number_input("Input desired spetrum number" , 
        #                                  step = 1, min_value = 1, max_value = len(filt_dfs[1])) 
        
        # Make a function thart shows one spectra at a time 
        # select one spectra from dataframe to pass into this fucntion
        
        #show_specta(Spectrum_select)
    

 
        
    


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




    
    
    

if __name__ == "__main__":
    main()
