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

    title = st.title("Internal ions")

    text_1 = st.markdown("desc")

    spectrum_file = st.file_uploader("Upload spectrum file:",
                                     type = ["raw"],
                                     help = "desc"
                                    )
      

    identifications_file = st.file_uploader("Upload identification file:",
                                            type = ["mzId"],
                                            help = "desc"
                                           )
    
    # RUN fraggannot?
    run_button = st.button("Analyze Fragments!", help = "Run fragannot")  
    
    if run_button:  
        json_data = run_fraggannot(spectrum_file, identifications_file)  
        pandas_df = json_pandas(json_data) 
        
        seq_length = st.number_input("Minimum sequence lenght")
        
        filt_df = filter_json_file(pandas_df, seq_length)
        

    if st.button("Make some plots!", help = "desc"):
        plot = st.plotly_chart(example_plot(), use_container_width = True) 
        
    
    
        

# side bar and main page loader
def main():

    about_str = \
    """
    description
    """

    st.set_page_config(page_title = "title",
                       page_icon = ":test_tube:",
                       layout = "centered",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/.../discussions",
                                     "Report a bug": "https://github.com/.../issues",
                                     "About": about_str}
                       )

    title = st.sidebar.title("title")

    logo = st.sidebar.image("img/logo.png", caption = "text")

    doc = st.sidebar.markdown("desc") 
    
    main_page()
    

        
    
    
    
def run_fraggannot(spectrum_file, identifications_file): 
    """
    Funtion that runs fragannot and outputs a json file.  
    
    
    
    """
    
    json_file = "file"

    return json_file

def json_pandas(json_file): 
    """ 
    
    """
    pd_df = "file"
    
    return pd_df


def filter_json_file(pd_df, int_frag_len): 
    """ 
    
    """
    filter_df = "file"
    
    return filter_df 



    
    
    

if __name__ == "__main__":
    main()
