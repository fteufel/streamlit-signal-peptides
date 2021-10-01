import streamlit as st
# To make things easier later, we're also importing numpy and pandas for
# working with sample data.
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from online_analysis_utils import region_features_from_server_output, make_one_dataframe
import altair as alt
import logomaker
#import plotly.express as px

st.title('Signal peptide explorer')


# load data
# @st.cache
# def load_data():
#     '''
#     Simple fn to load all archaea sequences.
#     '''
#     files_to_process = os.listdir('archaea_preds/')
#     print('Reading data...')

#     all_seqs = []
#     all_sizes = {} # for each proteome, also get total size (including negatives)
#     for f in files_to_process:
#         seqs = region_features_from_server_output(f'archaea_preds/{f}/')
#         seqs['Sequence'] = make_one_dataframe(f'archaea_preds/{f}/region_output.gff3', f'archaea_preds/{f}/output.gff3', f'archaea_preds/{f}/processed_entries.fasta')['Sequence']
#         seqs['Proteome_ID'] = f
        
#         df = pd.read_table('archaea_preds/'+ f +'/prediction_results.txt', skiprows=1)
#         all_sizes[f] = len(df)
        
#         all_seqs.append(seqs)
        
#     df = pd.concat(all_seqs)
#     return df, all_sizes



# Some plotting functions

def align_at_c_terminus(df):
    df_sp = df.loc[df['SP type'] == 'signal_peptide']
    max_len = int(df_sp['len_c'].max())
    aligned_seqs = []

    for idx, row in df_sp.iterrows():
        start = int(row['len_sp'] - row['len_c'])
        seq = row['Sequence'][start:row['len_sp']]
        #print(row['Sequence'][row['len_sp']:row['len_sp']+7])
        pad = '-' * int((max_len - len(seq)))
        aligned_seqs.append(pad + seq + row['Sequence'][row['len_sp']:row['len_sp']+7]) #add after padding

    return aligned_seqs


def make_logo(seqs):
    print('Plotting logo...')
    max_len = len(seqs[0]) -7 #because we add 7 after CS
    matrix = logomaker.alignment_to_matrix(seqs, characters_to_ignore='-', to_type='information')

    ww_logo = logomaker.Logo(matrix,
                            color_scheme='weblogo_protein',
                            vpad=.1,
                            width=.8)

    ww_logo.ax.set_xticks(list(range(0, max_len+7)))
    #plt.xlabels(list(range(-int(max_len)+1, 0, 1)) + list(range(1,8)))
    ww_logo.ax.set_xticklabels(list(range(-int(max_len)+1, 1, 1)) + list(range(1,8)))
    fig = plt.gcf()
    return fig


# Actual app.

id_name_df = pd.read_table('id_name_map.txt').set_index('Proteome_ID')

# df, all_sizes = load_data()
import json
with open('archaea_counts.json', 'r') as f:
    all_sizes = json.load(f)
df = pd.read_csv('archaea_data.csv', index_col=0)

option = st.sidebar.selectbox(
        'Choose a reference proteome',
        list(all_sizes.keys())
)

subset_df = df.loc[df['Proteome_ID']==option]
st.write('You selected:', f'{option}: {id_name_df.loc[option]["Species Name"]}')



st.header('Summary statistics')
st.text(f'Proteome size: {all_sizes[option]} sequences')
source =  pd.DataFrame(subset_df['SP type'].value_counts()).reset_index()
source =  pd.DataFrame(subset_df['SP type'].value_counts()).reset_index()
source = source.rename({'index': 'Type', 'SP type': 'Count'},axis=1)
source['Frequency (%)'] = source['Count'] / all_sizes[option] * 100
c = alt.Chart(source).mark_bar().encode(x = 'Type', y='Count', tooltip=['Frequency (%)']).interactive()
#c = alt.Chart(source).mark_bar().encode(y='SP type', x='index')
st.altair_chart(c, use_container_width=True)

# pd.Series({'Type''Count'=all_sizes[option]- source['Count'].sum(), 'Frequency (%):})
#fig = px.pie(df, values=source['Count'], names=source['Type'], title='Total proteome')
#st.plotly_chart(fig)


#st.write(subset_df)
st.header('Alignment at cleavage site')
seqs = align_at_c_terminus(subset_df)
fig = make_logo(seqs)
st.pyplot(fig)




