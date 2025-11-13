
import numpy as np
import scanpy as sc
import scdemon as sm
from scdemon.utils import recipe_full
from scdemon import plotting as pl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import scanpy as sc
import numpy as np
import diopy as dp
import os
import igraph as ig


data_py = dp.input.read_h5("Mal_obj_annotation.h5")
adata = data_py

prefixes = ['RTS','RTL','MT']
genes_to_remove = list(filter(lambda item: any(item.startswith(prefix) for prefix in prefixes),  adata.var_names))
print(genes_to_remove)

mask = ~adata.var_names.isin(genes_to_remove)
mask
adata = adata[:, mask]
adata.var_names_make_unique()

print(f"Remaining genes: {adata.shape[1]}")

recipe_full(adata, preprocess=True, annotate=True)

mod = sm.modules(adata, suffix='ESCC_example', k=100)
mod.setup() # Setup the object

print(pd.DataFrame(adata.obsm['X_umap'] , columns=['umap_1', 'umap_2'])  .head(10))

graph_id = 'base'
mod.make_graph(graph_id, resolution=2.5)

mlist = mod.get_modules(graph_id, print_modules=False)
mod.save_modules(graph_id)
scores = mod.graphs[graph_id].scores["leiden"]

np.save('scores.npy', scores)
scores_df = pd.DataFrame(scores)
scores_df.to_csv("scores_df.txt",sep='\t', index=False)

pl.plot_genes(mod, graph_id, attr="leiden", show_labels=True,ext='pdf', width=100)
pl.plot_genes(mod, graph_id, basis='umap', attr="leiden",ext='pdf', width=16)

pl.plot_umap_grid(mod, graph_id,plotname='plot_umap_grid',imgdir='./',width=3)

pl.plot_gene_logfc(mod, graph_id,attr="leiden",basis='umap',imgdir='./',width=16 , ext='pdf')


pl.plot_df_enrichment(mod, graph_id,df,ext='pdf')

pl.plot_heatmap_avgexpr(mod, graph_id,covariate_list=["sample.type", "Treatment"], cbar=True,attr="leiden",ext='pdf')

data_dict = mod.graphs[graph_id].__dict__ 

print(type(data_dict['mnames']))
print(data_dict['mnames'])
df_mnames = pd.DataFrame(data_dict['mnames'])
print(df_mnames.head())
df_mnames.to_csv('df_mnames.csv')

data_dict = mod.graphs[graph_id].__dict__

ig.write(data_dict['graph'], "graph.graphml")

graph = data_dict['graph']

edge_list = graph.get_edgelist()
vertex_names = graph.vs['name']

edges = []
for edge in edge_list:
  source_name = vertex_names[edge[0]]
target_name = vertex_names[edge[1]]
edges.append({'source': source_name, 'target': target_name})

edges_df = pd.DataFrame(edges)

edges_df.to_csv("edges.csv", index=False)
edges_df.to_csv("edges.txt",sep='\t', index=False)

print(edges_df.head())

gpres = sm.get_goterms(mod, graph_id)

for key, value in gpres.items():
  if isinstance(value, (pd.DataFrame, np.ndarray)):
  print(f"Key '{key}' has a non-serializable value: {type(value)}")

import os
target_path = os.path.join('', "goterms")

os.makedirs(target_path, exist_ok=True)
os.chdir(target_path)
for key, value in gpres.items():
  if isinstance(value, pd.DataFrame):
  filename = f"dataframe_{key}.csv"
value.to_csv(filename, index=False, encoding='utf-8')
print(f"Saved DataFrame '{key}' to '{filename}'")


import pandas as pd 
data_dict = mod.graphs[graph_id].__dict__
import numpy as np  

scores_array = data_dict['scores']
data_dict.keys()

print(f"Shape of 'scores': {scores_array['leiden'].shape}")
print(f"Data type of 'scores': {scores_array['leiden'].dtype}")

corr_array = data_dict['corr']   
corr_df = pd.DataFrame(corr_array)
corr_df.to_csv("corr_df.txt",sep='\t', index=False)


print("First few entries in 'scores':")
if scores_array['leiden'].ndim == 1:
  print(scores_array['leiden'][:5])
elif scores_array['leiden'].ndim == 2:
  print(scores_array['leiden'][:5, :5])
else:
  print("The array has more than 2 dimensions.")


def print_dict_info(d, level=0):
  """递归打印字典的结构、数据类型和内容。"""
if isinstance(d, dict):
  print(f"{'  '*level}Dictionary (level {level}):")
print(f"{'  '*level}Type: {type(d)}")

for key, value in d.items():
  print(f"{'  '*level}Key: {key}")
print_dict_info(value, level + 1)

elif isinstance(d, pd.DataFrame):
  print(f"{'  '*level}DataFrame (level {level}):")
print(f"{'  '*level}Type: {type(d)}")
print(f"{'  '*level}First two rows:")
print(d.head(2))

elif isinstance(d, list):
  print(f"{'  '*level}List (level {level}):")
print(f"{'  '*level}Type: {type(d)}")
print(f"{'  '*level}First few entries: {d[:2]}")

elif isinstance(d, np.ndarray):
  print(f"{'  '*level}Numpy Array (level {level}):")
print(f"{'  '*level}Type: {type(d)}")
print(f"{'  '*level}Shape: {d.shape}")
print(f"{'  '*level}First few entries: {d.flatten()[:2]}")

else:
  print(f"{'  '*level}Value (level {level}): {d} (Type: {type(d)})") 


print_dict_info(gpres)



