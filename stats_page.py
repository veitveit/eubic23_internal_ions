import json
import pandas as pd
import plotly.io as io
import numpy as np

# 
import plotly.graph_objects as go
import plotly.figure_factory as ff




def example_plot():
    x1 = np.random.randn(200) - 2
    x2 = np.random.randn(200)
    x3 = np.random.randn(200) + 2
    hist_data = [x1, x2, x3]
    group_labels = ['Group 1', 'Group 2', 'Group 3']
    fig = ff.create_distplot(hist_data, group_labels, bin_size = [.1, .25, .5])
    
    return fig



#####  ------DATA INPUT-------
#json_file = "C:/Users/0kiko/Downloads/Human-Protein-Training_Trypsin.json"
#with open(json_file, "r", encoding = "utf-8") as f:
#    data = json.load(f)
#    f.close()

#dataframes = json_to_dataframes("C:/Users/0kiko/Downloads/Human-Protein-Training_Trypsin.json")

#fragments_dataframe = dataframes[0]
#spectra_dataframe = dataframes[1]

##### filtering

#### Visualizations ----DATA PROCESSING?----
def common_type_gen(fragments_dataframe):
    common_type = fragments_dataframe['frag_type1'].astype(str).str.cat(fragments_dataframe['frag_type2'], sep='-')
    common_type = common_type.replace("n", "not annotated", regex=True)
    common_type = common_type.replace("t-","",regex=True)
    common_type = common_type.replace("-t","",regex=True)
    #counts = common_type.value_counts()
    fragments_dataframe['frag_types'] = common_type
    return common_type
    
    
    
    
# Frequency of ion types 
def common_type_hist(fragments_dataframe):
    
    common_type = common_type_gen(fragments_dataframe)
    
    fig = go.Figure([go.Histogram(x=common_type)])
    return fig



    
def common_type_pie(fragments_dataframe):
    
    common_type = common_type_gen(fragments_dataframe)
    
    counts = common_type.value_counts()
    fig = go.Figure([go.Pie(labels=counts.keys(), values=counts)])
    return fig




# intensity distribution of different ions
# density or probability
def log_ion_intens_dist(fragments_dataframe):
    histnorm = "probability"
    types = fragments_dataframe["frag_types"].unique() 
    #print(types)
    histograms = list()
    for t in types:
        histograms.append(go.Histogram(x=np.log(fragments_dataframe[fragments_dataframe.frag_types == t].frag_intensity),histnorm=histnorm, name=t, nbinsx=50))
        fig = go.Figure(histograms)
        fig.update_layout(
            barmode='group',
            title="Histograms of logarithmic intensities per ion type",
            xaxis_title="log2(intensity)",
            yaxis_title=histnorm)
    return fig



# relative intensity to total intensity distribution of different ions
# density or probability
def rel_ion_intens_perc(fragments_dataframe):
    histnorm = "probability"
    types = fragments_dataframe["frag_types"].unique() 
    #print(types)
    histograms = list() 
    for t in types:
        histograms.append(go.Histogram(x=fragments_dataframe[fragments_dataframe.frag_types == t].perc_of_total_intensity,histnorm=histnorm, name=t, nbinsx=50))
        fig = go.Figure(histograms)
        fig.update_layout(
            barmode='group',
            title="Histograms of relative intensities per ion type",
            xaxis_title="intensity",
            yaxis_title=histnorm)

    return fig



# relative intensity to total intensity distribution of different ions
# density or probability
def rel_ion_intens_prop(fragments_dataframe):
    histnorm = "probability"
    types = fragments_dataframe["frag_types"].unique() 
    #print(types)
    histograms = list()
    for t in types:
        histograms.append(go.Histogram(x=fragments_dataframe[fragments_dataframe.frag_types == t].prop_intensity_to_base_peak,histnorm=histnorm, name=t, nbinsx=50))
        fig = go.Figure(histograms)
        fig.update_layout(
            barmode='group',
            title="Histograms of relative intensities to base peak per ion type",
            xaxis_title="Percentage per base peak",
            yaxis_title=histnorm)
        #fig.show()
    return fig


# mz distribution of different ion types
# density or probability
def mz_dist_ion_type(fragments_dataframe):
    histnorm = "probability"
    types = fragments_dataframe["frag_types"].unique() 
    #print(types)
    histograms = list()
    for t in types:
        histograms.append(go.Histogram(x=fragments_dataframe[fragments_dataframe.frag_types == t].frag_mz, histnorm=histnorm, name=t, nbinsx=50))
        fig = go.Figure(histograms)
        fig.update_layout(
            barmode='group',
            title="Histograms of mz values per ion type",
            xaxis_title="mz",
            yaxis_title=histnorm)
        #fig.show()
    return fig



### Percentages of different ion types per spectrum
def per_spec_ion_type(spectra_dataframe):
    histnorm = "probability"
    types = ["internal","terminal","other"]
    print(types)
    histograms = list()
    for t in types:
        histograms.append(go.Histogram(x=spectra_dataframe["perc_" + t], histnorm=histnorm, name=t, nbinsx=50))
        fig = go.Figure(histograms)
        fig.update_layout(
            barmode='group',
            title="Histograms of percentages of ion type per spectrum",
            xaxis_title="Percentage",
            yaxis_title=histnorm)
        #fig.show()
    return fig




### Same for intensities
def per_spec_ion_intens(spectra_dataframe):
    histnorm = "probability"
    types = ["internal","terminal","other"]
    histograms = list()
    for t in types:
        histograms.append(go.Histogram(x=spectra_dataframe["total_int_" + t], histnorm=histnorm, name=t, nbinsx=50))
        fig = go.Figure(histograms)
        fig.update_layout(
            barmode='group',
            title="Histograms of percentages of ion type per spectrum",
            xaxis_title="Percentage",
            yaxis_title=histnorm)
        #fig.show()
    return fig




### Ridgelines of logarithmic intensities per ion type
def log_ion_intens_ridge(fragments_dataframe):
    histnorm = "probability"
    types = fragments_dataframe["frag_types"].unique()
    fig = go.Figure()
    for t in types:
        fig.add_trace(go.Violin(x=np.log(fragments_dataframe[fragments_dataframe.frag_types == t].frag_intensity),name=t))
        fig.update_traces(orientation='h', side='positive', width=3, points=False)
        fig.update_layout(
            barmode='group',
            title="Ridgelines of logarithmic intensities per ion type",
            xaxis_title="log2(intensity)",
            yaxis_title=histnorm)
    return fig




### Ridgelines of relative intensities per ion type
def rel_ion_intens_ridge(fragments_dataframe):
    fig = go.Figure()
    histnorm = "probability"
    types = fragments_dataframe["frag_types"].unique()
    for t in types:
        fig.add_trace(go.Violin(x=fragments_dataframe[fragments_dataframe.frag_types == t].perc_of_total_intensity,name=t))
        fig.update_traces(orientation='h', side='positive', width=3, points=False)
        fig.update_layout(
            barmode='group',
            title="Ridgelines of relative intensities per ion type",
            xaxis_title="intensity",
            yaxis_title=histnorm)
    return fig




## Ridgelines of relative intensities to base peak per ion type
def rel_ion_intens_prop_ridge(fragments_dataframe):
    histnorm = "probability"
    types = fragments_dataframe["frag_types"].unique()
    fig = go.Figure()
    for t in types:
        fig.add_trace(go.Violin(x=fragments_dataframe[fragments_dataframe.frag_types == t].prop_intensity_to_base_peak,name=t))
        fig.update_traces(orientation='h', side='positive', width=3, points=False)
        fig.update_layout(
            barmode='group',
            title="Ridgelines of relative intensities to base peak per ion type",
            xaxis_title="intensity",
            yaxis_title=histnorm)
    return fig






## Distribution of total intensity of single amino acids
# Filter for all single amino acid fragments
#fragments_dataframe[fragments_dataframe['frag_seq'].str.len() == 1 and fragments_dataframe.modification.empty()]