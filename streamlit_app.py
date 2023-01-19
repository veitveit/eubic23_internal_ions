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

from json_to_dataframes import json_to_dataframes, json_to_dataframes2
from fragannot_call import fragannot_call
from stats_page import *  # plotting functions
from spectra_page import *  # spectra view
from filtering_files import *  # filter files functions




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
        formula = st.latex(r"""M(fragment)= \frac{M(peptide)+\Delta M(IonCap_{start})-\Delta M(IonCap_{end})+M(H)\times charge -M(formular) } {charge}""")

    with data_import_tab:

        title = st.title("Data Import")

        spectrum_file = st.file_uploader("Upload spectrum file:",
                                         type = ["raw", "mgf"],
                                         help = "desc"
                                        )

        identifications_file = st.file_uploader("Upload identification file:",
                                                type = ["mzid"],
                                                help = "desc"
                                               )

        tolerance = st.number_input("Tolerance in Da",
                                    value = 0.02)
        fions_text = st.markdown("Select which ion types of your MS run:")

        # List of booleans

        fragannot_ions_col1, fragannot_ions_col2 = st.columns(2)

        with fragannot_ions_col2:

            fions_checkbox_nterm = st.markdown("N-terminal ions")

            fA_ion = st.checkbox("A ions", key = "fA_ion")
            fB_ion = st.checkbox("B ions", key = "fB_ion")
            fC_ion = st.checkbox("C ions", key = "fC_ion")
            fCdot_ion = st.checkbox("Cdot ions", key = "fCdot_ion")
            fCm1_ion = st.checkbox("C-1 ions", key = "fCm1_ion")
            fCp1_ion = st.checkbox("C+1 ions", key = "fCp1_ion")

        with fragannot_ions_col1:

            fions_checkbox_cterm = st.markdown("C-terminal ions")

            fX_ion = st.checkbox("X ions", key = "fX_ion")
            fY_ion = st.checkbox("Y ions", key = "fY_ion")
            #fZ_ion = st.checkbox("Z ions", key = "fZ_ion")
            fZdot_ion = st.checkbox("Zdot ions", key = "fZdot_ion")
            fZp1_ion = st.checkbox("Z+1 ions", key = "fZp1_ion")
            fZp2_ion = st.checkbox("Z+2 ions", key = "fZp2_ion")
            fZp3_ion = st.checkbox("Z+3 ions", key = "fZp3_ion")

        #fragannot_ion_selection = [fN_ion, fA_ion, fB_ion, fC_ion, fCdot_ion, fCm1_ion, fCp1_ion, fX_ion, fY_ion, fZ_ion, fZdot_ion, fZp1_ion, fZp2_ion, fZp3_ion]
        fragannot_ion_selection = [fA_ion, fB_ion, fC_ion, fCdot_ion, fCm1_ion, fCp1_ion, fX_ion, fY_ion, fZdot_ion, fZp1_ion, fZp2_ion, fZp3_ion]
        fragannot_ion_names = ["a", "b", "c", "cdot", "c-1", "c+1", "x", "y", "zdot", "z+1", "z+2", "z+3"]
        fragannot_call_ion_selection = []
        for i, sel in enumerate(fragannot_ion_selection):
            if sel:
                fragannot_call_ion_selection.append(fragannot_ion_names[i])

        charges_str = st.text_input("Charges to consider [comma delimited]",
                                    value = "-1,+1")
        charges = charges_str.split(",")
        losses_str = st.text_input("Neutral losses to consider [comma delimited]",
                                   value = "H2O")
        losses = losses_str.split(",")

        data_text = st.markdown("Upload and filter  json files with internal fragment ion matches from fragannot here. https://github.com/arthur-grimaud/fragannot")

        test_FILE = st.file_uploader("Upload json file:",
                                                type = ["json"],
                                                help = "Process you mzid and mgf file to fins internal fragment ions using fragannot on: https://github.com/arthur-grimaud/fragannot"
                                               )

        data_text = st.markdown("Select filtering options below")

        start_seq_length, end_seq_length = st.select_slider("Peptide sequence lenght",
                                            options= range(0,5001),
                                            value = (0,5000))

        start_frag_len, end_frag_len = st.select_slider("Fragment ion sequence lenght",
                                            options= range(0,1001),
                                            value = (0,1000))

        start_mz, end_mz = st.select_slider("Select m/z range for internal fragment",
                                            options= range(0,100001),
                                            value = (0,100000))

        start_int, end_int = st.select_slider("Select intensity range for internal fragment ",
                                              options= range(0,1000001),
                                              value =(0,1000000))

        #seq_length = st.number_input("Identification score range")
        data_text = st.markdown("Select which ions to analyse")

        N_ion = st.checkbox("Non-annotated ions", value = True)

        # List of booleans

        check_col1, check_col2 = st.columns(2)

        with check_col2:

            ions_checkbox_nterm = st.markdown("N-terminal ions")

            A_ion = st.checkbox("A ions", value = True)
            B_ion = st.checkbox("B ions", value = True)
            C_ion = st.checkbox("C ions", value = True)
            Cdot_ion = st.checkbox("Cdot ions", value = True)
            Cm1_ion = st.checkbox("C-1 ions", value = True)
            Cp1_ion = st.checkbox("C+1 ions", value = True)

        with check_col1:

            ions_checkbox_cterm = st.markdown("C-terminal ions")

            X_ion = st.checkbox("X ions", value = True)
            Y_ion = st.checkbox("Y ions", value = True)
            Z_ion = st.checkbox("Z ions", value = True)
            Zdot_ion = st.checkbox("Zdot ions")
            Zp1_ion = st.checkbox("Z+1 ions", value = True)
            Zp2_ion = st.checkbox("Z+2 ions", value = True)
            Zp3_ion = st.checkbox("Z+3 ions", value = True)

        ion_filter_param = [N_ion, A_ion, B_ion, C_ion, Cdot_ion, Cm1_ion, Cp1_ion, X_ion, Y_ion, Z_ion, Zdot_ion, Zp1_ion, Zp2_ion, Zp3_ion]

        l1, l2, center_button, r1, r2 = st.columns(5)

        with center_button:
            run_analysis = st.button("Analyze files!")

        if test_FILE is not None or (spectrum_file is not None and identifications_file is not None):

            if run_analysis:
                if test_FILE is not None:
                    dataframes = json_to_dataframes(test_FILE, is_file = True)
                else:
                    fragannot_result = fragannot_call(spectrum_file,
                                                      identifications_file,
                                                      float(tolerance),
                                                      fragannot_call_ion_selection,
                                                      charges,
                                                      losses)
                    dataframes = json_to_dataframes2(fragannot_result)




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

                    col111, col222 = st.columns(2)

                    with col111:
                        title = st.header("Ion type per spectra")
                        plot = st.plotly_chart(per_spec_ion_type(filt_dfs[1]), use_container_width = True)

                    with col222:
                        title = st.header("Log Intensities")
                        plot = st.plotly_chart(per_spec_ion_intens(filt_dfs[1]), use_container_width = True)

                    # Logo view
                    title = st.title("Logo view of internal fragments")
                    #plot = st.plotly_chart(, use_container_width = True)

                    plot = st.pyplot(logo_of_fraction(filt_dfs[1])) # , clear_figure=None, **kwargs


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

    logo = st.sidebar.image("img/image(1).gif") # , caption = "Logo"

    doc = st.sidebar.markdown("desc")


    contact_str = "**Contact:** [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk), [Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk), [Kristian Boje Nielsen](krisn16@student.sdu.dk), [Louise Buur](louise.buur@fh-hagenberg.at), [Micha Birklbauer](mailto:micha.birklbauer@gmail.com), [Mohieddin Jafari](mohieddin.jafari@helsinki.fi), [Veit S](veits@bmb.sdu.dk), [Vladimir Gorshkov](homer2k@gmail.com), [Zoltan Udvardy](zoltan.udvardy.ipbs@gmail.com) "
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/michabirklbauer/piaweb/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)


    main_page()









if __name__ == "__main__":
    main()
