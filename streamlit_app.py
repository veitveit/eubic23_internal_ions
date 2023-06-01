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

from utils.stats_plotting import *  # plotting functions
from utils.spectra_plotting import *  # spectra view
from utils.filter_dataframes import filter_dataframes  # filter files functions

# main page content
def main_page():

    title = st.title("Fragment Explorer")

    sub_title1 = st.header("About Fragannot and Fragment Explorer")

    general_description = \
    """
    We created a workflow for exploring internal ions from raw spectra and identified spectra given in various formats. A comprehensive nomenclature for internal ions allows their precise definition and the calculation of their masses. The nomenclature was implemented into a new tool for annotating fragment ions that we titled **"Fragannot"**. This tool outputs fragment annotations as a *.json* or *.csv* files. The *.csv* files can then read with the tool - **"Fragment Explorer"** - to create fragment centric and spectrum centric statistics as well as multiple interactive visualizations.

    Fragannot can be downloaded for free from [github.com/arthur-grimaud/fragannot](https://github.com/arthur-grimaud/fragannot).
    """
    text_1 = st.markdown(general_description)

    # Decripbe project, goal and members
    text_2 = st.markdown("We are calculating the mass of the fragments as described in the following:")

    # Add forumula for nomenclatur and mass
    formula = st.latex(r"""M(fragment)= \frac{M(peptide)+\Delta M(IonCap_{start})-\Delta M(IonCap_{end})+M(H)\times charge -M(formular) } {charge}""")

    text_3 = st.markdown("Fragannot results can be analyzed below:")

    sub_title2 = st.header("Data Import")

    fragments = st.file_uploader("Upload the fragment-centric .csv output of Fragannot:",
                                 type = ["csv"],
                                 help = "Upload the fragment-centric .csv output of Fragannot."
                                )

    spectra = st.file_uploader("Upload the spectrum-centric .csv output of Fragannot:",
                               type = ["csv"],
                               help = "Upload a identification file that contains PSMs of the spectrum file in mzID format."
                              )

    text_4 = st.markdown("Fragannot results can be further filtered with the filtering options below (optional).")

    with st.expander("Expand filtering options:"):

        text_5 = st.markdown("**Select filtering options below:**")

        start_seq_length, end_seq_length = st.select_slider("Peptide sequence lenght:",
                                                            options = range(0, 5001),
                                                            value = (0, 5000),
                                                            help = "For 'Spectrum-centric Statistics' only spectra with peptides withing range are considered.")

        start_frag_len, end_frag_len = st.select_slider("Fragment ion sequence lenght:",
                                                        options = range(0, 1001),
                                                        value = (0, 1000),
                                                        help = "For 'Fragment-centric Statistics' only fragments within in range are considered.")

        start_mz, end_mz = st.select_slider("Select m/z range for internal fragments:",
                                            options = range(0, 100001),
                                            value = (0, 100000),
                                            help = "For 'Fragment-centric Statistics' only fragments within in range are considered.")

        start_int, end_int = st.select_slider("Select intensity range for internal fragments:",
                                              options = range(0, 1000001),
                                              value =(0, 1000000),
                                              help = "For 'Fragment-centric Statistics' only fragments within in range are considered.")

        data_text = st.markdown("Select which ions to analyse.")

        N_ion = st.checkbox("Non-annotated ions", value = True)

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
            Zdot_ion = st.checkbox("Zdot ions", value = True)
            Zp1_ion = st.checkbox("Z+1 ions", value = True)
            Zp2_ion = st.checkbox("Z+2 ions", value = True)
            Zp3_ion = st.checkbox("Z+3 ions", value = True)

        ion_filter_param = [N_ion, A_ion, B_ion, C_ion, Cdot_ion, Cm1_ion, Cp1_ion, X_ion, Y_ion, Z_ion, Zdot_ion, Zp1_ion, Zp2_ion, Zp3_ion]

    l1, l2, center_button, r1, r2 = st.columns(5)

    with center_button:
        run_analysis = st.button("Analyze and filter files!")

    if run_analysis and fragments is not None and spectra is not None:
        st.session_state["fragments_df"] = pd.read_csv(fragments)
        st.session_state["spectra_df"] = pd.read_csv(spectra)
        st.session_state["filtered_dfs"] = filter_dataframes([st.session_state["fragments_df"], st.session_state["spectra_df"]],
                                                             start_seq_length,
                                                             end_seq_length,
                                                             start_frag_len,
                                                             end_frag_len,
                                                             start_mz,
                                                             end_mz,
                                                             start_int,
                                                             end_int,
                                                             ion_filter_param)

        data_text = st.success("Finished data processing!")

    if "filtered_dfs" in st.session_state:
        filt_dfs = st.session_state["filtered_dfs"]

        st.table(filt_dfs[0].head(10))

        sub_title3 = st.header("Fragment-centric Statistics")

        plot1_col1, plot2_col2 = st.columns(2)

        with plot1_col1:
            plot1_title = st.header("Histogram of ion types:")
            plot1 = st.plotly_chart(common_type_hist(filt_dfs[0]), use_container_width = True)

        with plot2_col2:
            plot2_title = st.header("Pie chart of ion types:")
            plot2 = st.plotly_chart(common_type_pie(filt_dfs[0]), use_container_width = True)

        plot3_title = st.header("Histogram of m/z per ion type:")
        plot3 = st.plotly_chart(mz_dist_ion_type(filt_dfs[0]), use_container_width = True)

        plot4_col1, plot5_col2 = st.columns(2)

        with plot4_col1:
            plot4_title = st.header("Log Intensities Distribution:")
            plot4 = st.plotly_chart(rel_ion_intens_perc(filt_dfs[0]), use_container_width = True)

        with plot5_col2:
            # relative intensity to total intensity distribution of different ions
            plot5_title = st.header("Relative Log Intensities:")
            plot5 = st.plotly_chart(rel_ion_intens_ridge(filt_dfs[0]), use_container_width = True)

        sub_title4 = st.header("Spectrum-centric Statistics")

        plot6_col1, plot7_col2 = st.columns(2)

        with plot6_col1:
            plot6_title = st.header("Ion type per spectrum:")
            plot6 = st.plotly_chart(per_spec_ion_type(filt_dfs[1]), use_container_width = True)

        with plot7_col2:
            plot7_title = st.header("Log Intensities")
            plot7 = st.plotly_chart(per_spec_ion_intens(filt_dfs[1]), use_container_width = True)

        plot8_title = st.header("Logo view of internal fragments:")

        plot8_text_1 = st.markdown("Select number of top spectra with the highest number of internal ions to inlcude in logo. Choose 0 for making a logo of all spectra.")
        plot8_topn = st.number_input("Choose number of spectra in logo:", min_value = 0, max_value = 3, step = 1)
        p8_min_length_logo, p8_max_length_logo = st.select_slider("Peptide sequence lenght:",
                                                            options = range(0,21),
                                                            value = (0,20))
        plot8 = st.pyplot(logo_of_fraction(filt_dfs[1], filt_dfs[0], plot8_topn, p8_max_length_logo, p8_min_length_logo))

# side bar and main page loader
def main():

    about_str = \
    """
    Fragment Explorer is a small program to explore fragment ions in mass spectra.
    """

    st.set_page_config(page_title = "Internal ions",
                       page_icon = ":test_tube:",
                       layout = "wide",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/veitveit/eubic23_internal_ions/discussions",
                                     "Report a bug": "https://github.com/veitveit/eubic23_internal_ions/issues",
                                     "About": about_str}
                       )

    title = st.sidebar.title("Fragment Explorer - EuBIC 2023 Hackathon")

    logo = st.sidebar.image("img/logo.png")

    doc = st.sidebar.markdown(about_str)

    contact_str = "**Contact:** [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk), [Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk), [Kristian Boje Nielsen](krisn16@student.sdu.dk), [Louise Buur](louise.buur@fh-hagenberg.at), [Micha Birklbauer](mailto:micha.birklbauer@gmail.com), [Mohieddin Jafari](mohieddin.jafari@helsinki.fi), [Veit Schwämmle](veits@bmb.sdu.dk), [Vladimir Gorshkov](homer2k@gmail.com), [Zoltan Udvardy](zoltan.udvardy.ipbs@gmail.com) "
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/veitveit/eubic23_internal_ions/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)

    main_page()

if __name__ == "__main__":
    main()
