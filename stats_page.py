import json
import pandas as pd
import plotly.io as io
import plotly.graph_objects as go
import numpy as np



def example_plot():
    x1 = np.random.randn(200) - 2
    x2 = np.random.randn(200)
    x3 = np.random.randn(200) + 2
    hist_data = [x1, x2, x3]
    group_labels = ['Group 1', 'Group 2', 'Group 3']
    fig = ff.create_distplot(hist_data, group_labels, bin_size = [.1, .25, .5])
    return fig



#io.renderers.default="browser"




### Summarize fragment ions over all spectra
def all_ion_sum(filename):
    f = open(filename)
    data = json.load(f)
    all_ions = pd.DataFrame()
    for i in data:
        all_ions = pd.concat([all_ions, pd.DataFrame.from_dict(i['annotation'])])
        #all_ions.info()
        #all_ions
### separate ion types, ..
    all_ions[['frag_types', 'frag_rest']] = all_ions.theoretical_code.str.split("@", expand = True)
    all_ions[['frag_type_1', 'frag_type_2']] = all_ions.frag_types.str.split(":", expand = True)
    all_ions[['frag_length', 'frag_rest']] = all_ions.frag_rest.str.split("(", expand = True)
    all_ions[['position_frag_type_1', 'position_frag_type_2']] = all_ions.frag_length.str.split(":", expand = True)
    all_ions['frag_length'] =  all_ions.position_frag_type_2.astype(float) - all_ions.position_frag_type_1.astype(float) + 1
    return all_ions




#### Visualizations
def common_type(all_ions):
    common_type = all_ions['frag_type_1'].astype(str).str.cat(all_ions['frag_type_2'], sep='-')
    common_type = common_type.fillna("not annotated")
    common_type = common_type.replace("t-","",regex=True)
    common_type = common_type.replace("-t","",regex=True)

    all_ions['frag_types'] = common_type
    return common_type




# Frequency of ion types
def common_type_hist(common_type):
    fig = go.Figure([go.Histogram(x=common_type)])
    fig.show()
    return fig




def common_type_pie(common_type):
    counts = common_type.value_counts()
    fig = go.Figure([go.Pie(labels=counts.keys(), values=counts)])
    fig.show()
    return fig



### intensity distribution of different ions
### density or probability
def intens_dist_ion(all_ions):
    histnorm = "probability"
    types = all_ions["frag_types"].unique() 
    histograms = []
    for t in types:
        histograms = histograms.append(go.Histogram(x=np.log(all_ions[all_ions.frag_type_1 == t].intensity),histnorm=histnorm))
        fig = go.Figure([histograms])
        #fig.show()
        return fig
    
    

### distribution of first ions
def first_ion_dist(all_ions):
    fig = go.Figure([go.Histogram(x=(all_ions[all_ions.frag_type_1 == "b"].mz),histnorm='probability'),
                     go.Histogram(x=(all_ions[all_ions.frag_type_1 == "t"].mz),histnorm='probability')])
    #fig.show()
    return fig
    

