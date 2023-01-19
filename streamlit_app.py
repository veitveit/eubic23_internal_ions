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

    main_tab, data_import_tab, stats_tab = st.tabs(["Main page",
                                                    "Data import",
                                                    "Statistical analysis"]) # "Spectra analysis",  spectra_tab,

    with main_tab:
        title = st.title("Fragment Exlporer")

        general_description = \
        """
        We created a workflow for exploring internal ions from raw spectra and identified spectra given in various formats. A comprehensive nomenclature for internal ions allows their precise definition and the calculation of their masses. The nomenclature was implemented into a new tool for annotating fragment ions entitled "fragannot". This tool output fragment annotation as a json file. This json file is then read by "fragment explorer" to create fragment centric and spectrum centric statistics as well as multiple interactive visualizations.
        """
        text_1 = st.markdown(general_description)

        # Decripbe project, goal and members
        text_1 = st.markdown("We are calculating the mass of the fragments as described in the following:")

        # Add forumula for nomenclatur and mass
        formula = st.latex(r"""M(fragment)= \frac{M(peptide)+\Delta M(IonCap_{start})-\Delta M(IonCap_{end})+M(H)\times charge -M(formular) } {charge}""")

    with data_import_tab:

        title = st.title("Data Import")

        spectrum_file = st.file_uploader("Upload spectrum file:",
                                         type = ["raw", "mgf"],
                                         help = "Upload a spectrum file to be analyzed in .raw of .mgf format."
                                        )

        identifications_file = st.file_uploader("Upload identification file:",
                                                type = ["mzid"],
                                                help = "Upload a identification file that contains PSMs of the spectrum file in mzID format."
                                               )

        tolerance = st.number_input("Tolerance in Da:",
                                    value = 0.02,
                                    help = "Fragment mass tolerance in Dalton.")
        fions_text = st.markdown("Select which ion types your applied fragmentation method produced:")

        # List of booleans

        fragannot_ions_col1, fragannot_ions_col2 = st.columns(2)

        with fragannot_ions_col2:

            fions_checkbox_nterm = st.markdown("**N-terminal ions:**")

            fA_ion = st.checkbox("A ions", key = "fA_ion")
            fB_ion = st.checkbox("B ions", value = True, key = "fB_ion")
            fC_ion = st.checkbox("C ions", key = "fC_ion")
            fCdot_ion = st.checkbox("Cdot ions", key = "fCdot_ion")
            fCm1_ion = st.checkbox("C-1 ions", key = "fCm1_ion")
            fCp1_ion = st.checkbox("C+1 ions", key = "fCp1_ion")

        with fragannot_ions_col1:

            fions_checkbox_cterm = st.markdown("**C-terminal ions:**")

            fX_ion = st.checkbox("X ions", key = "fX_ion")
            fY_ion = st.checkbox("Y ions", value = True, key = "fY_ion")
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

        charges_str = st.text_input("Charges to consider [comma delimited]:",
                                    value = "-1,+1",
                                    help = "The charges to consider for fragment ions. Multiple entries should be delimited by commas!")
        charges = charges_str.split(",")
        losses_str = st.text_input("Neutral losses to consider [comma delimited]",
                                   value = "H2O",
                                   help = "Neutral losses to consider for fragment ions. Multiple entries should be delimited by commas!")
        losses = losses_str.split(",")

        data_text = st.markdown("Upload and filter .json files with fragment ion annotations from fragannot here. https://github.com/arthur-grimaud/fragannot")

        json_file = st.file_uploader("Upload json file:",
                                     type = ["json"],
                                     help = "Process your .mzid and .mgf file to annotate fragment ions using fragannot on: https://github.com/arthur-grimaud/fragannot"
                                    )

        data_text = st.markdown("**Select filtering options below:**")

        start_seq_length, end_seq_length = st.select_slider("Peptide sequence lenght:",
                                            options = range(0,5001),
                                            value = (0,5000),
                                            help = "For 'Spectrum-centric Statistics' only spectra with peptides withing range are considered.")

        start_frag_len, end_frag_len = st.select_slider("Fragment ion sequence lenght:",
                                            options = range(0,1001),
                                            value = (0,1000),
                                            help = "For 'Fragment-centric Statistics' only fragments within in range are considered.")

        start_mz, end_mz = st.select_slider("Select m/z range for internal fragments:",
                                            options = range(0,100001),
                                            value = (0,100000),
                                            help = "For 'Fragment-centric Statistics' only fragments within in range are considered.")

        start_int, end_int = st.select_slider("Select intensity range for internal fragments:",
                                              options = range(0,1000001),
                                              value =(0,1000000),
                                              help = "For 'Fragment-centric Statistics' only fragments within in range are considered.")

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

        if json_file is not None or (spectrum_file is not None and identifications_file is not None):

            if run_analysis:
                if json_file is not None:
                    dataframes = json_to_dataframes(json_file, is_file = True)
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
                data_text = st.markdown("Finished data processing!")
                st.table(filt_dfs[0].head(10))

                with stats_tab:

                    title = st.title("Fragment-centric Statistics")

                    col1, col2 = st.columns(2)

                    with col1:
                        title = st.header("Histogram of ion types:")
                        plot = st.plotly_chart(common_type_hist(filt_dfs[0]), use_container_width = True)

                    with col2:
                        title = st.header("Pie chart of ion types:")
                        plot = st.plotly_chart(common_type_pie(filt_dfs[0]), use_container_width = True)

                    title = st.header("Histogram of m/z per ion type:")
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

                    per_spectrum_title = st.title("Spectrum-centric Statistics")

                    with col111:
                        title = st.header("Ion type per spectrum")
                        plot = st.plotly_chart(per_spec_ion_type(filt_dfs[1]), use_container_width = True)

                    with col222:
                        title = st.header("Log Intensities")
                        plot = st.plotly_chart(per_spec_ion_intens(filt_dfs[1]), use_container_width = True)

                    # Logo view
                    title = st.header("Logo view of internal fragments")
                    #plot = st.plotly_chart(, use_container_width = True)

                    #text_1 = st.markdown("Select number of top spectra with the highest number of internal ions to inlcude in logo. Choose number 0 for making a logo of all spectra.")
                    #number_topn = st.number_input("Choose number of spectra in logo:", min_value = 0, max_value = 3, step = 1)
                    #min_length_logo, max_length_logo = st.select_slider("Peptide sequence lenght:",
                    #                                                    options = range(0,21),
                    #                                                    value = (0,20))
                    plot = st.pyplot(logo_of_fraction(filt_dfs[1], filt_dfs[0]))#, number_topn, max_length_logo, min_length_logo))

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
    Fragment Explorer is a small program to explore fragment ions in mass spectra.
    """

    st.set_page_config(page_title = "Internal ions",
                       page_icon = "random", # ":test_tube:"
                       layout = "wide",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/veitveit/eubic23_internal_ions/discussions",
                                     "Report a bug": "https://github.com/veitveit/eubic23_internal_ions/issues",
                                     "About": about_str},  #

                       )

    title = st.sidebar.title("Fragment Explorer - EuBIC 2023 Hackathon")

    logo = st.sidebar.image("img/image(1).gif") # , caption = "Logo"

    doc = st.sidebar.markdown(about_str)

    contact_str = "**Contact:** [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk), [Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk), [Kristian Boje Nielsen](krisn16@student.sdu.dk), [Louise Buur](louise.buur@fh-hagenberg.at), [Micha Birklbauer](mailto:micha.birklbauer@gmail.com), [Mohieddin Jafari](mohieddin.jafari@helsinki.fi), [Veit Schwämmle](veits@bmb.sdu.dk), [Vladimir Gorshkov](homer2k@gmail.com), [Zoltan Udvardy](zoltan.udvardy.ipbs@gmail.com) "
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/veitveit/eubic23_internal_ions/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)

    main_page()

if __name__ == "__main__":
    main()
