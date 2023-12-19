# MALDIpy: Single-cell analysis of MALDI-MS imaging mass spectrometry data
`MALDIpy` is a Python package for analyzing MALDI-MS imaging mass spectrometry data at the single-cell level (each “cell” refers to a 10-µm pixel metabolome).

Its function includes: metabolite feature visualization, Scanpy-based single-cell analysis with UMAP clustering, projection of cluster annotations, integrative multi-sample visualization and analysis.

<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/function.png" width="600"> <br>

# 1. Installation

To install `MALDIpy`, run it in a terminal

```
pip install MALDIpy
```

# 2. Tutorials

(1) Analyze one single dataset with high efficiency -- all analysis can be done in 5-10 minutes on a labtop:
https://github.com/TheHumphreysLab/MALDIpy/blob/main/vignette/kidney_cortex.ipynb

(2) Integrative analysis of multiple MALDI-MS datasets:
https://github.com/TheHumphreysLab/MALDIpy/blob/main/vignette/kidney_integration.ipynb
<br>

# 3. Key functions

(1) Create a MALDIpy object with `msi_data`.

```
maldi_obj = msi_data(raw_file, scale=10)
```

(2) Visualize any metabolite of interest in a MALDIpy object with `plt`.

```
maldi_obj.plt(mz=741.530654593237, figsize = (6,5), smooth=False, pos = 'lower left', remove_hs = True, cmap = "magma_r")
```
<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/plt1.png" width="300"> <br>

(3) MALDIpy to AnnData conversion with `to_adata`.

```
adata = maldi_obj.to_adata(add_meta=True, csv_file=raw_file)
```

(4) Efficient single-cell quality control, dimension reduction and clustering with `single_cell`, including `single_cell.maldifilter`, `single_cell.maldi_norm` and `single_cell.maldi_clustering`.<br>
`single_cell.maldi_clustering` includes Harmony-based batch effect correction when processing multiple tissue sections.<br>
<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/plt2.png" width="300"><img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/plt4.png" width="450"> <br>

(5) Project single-cell cluster annotation onto the tissue section with `projection`.
```
MALDIpy.projection.umap_projection(adata, file_name=raw_file,pltcmap=adata.uns['leiden_colors'],
                                   figtitle='Leiden Cluster Projection',figdpi=150, fig_size=(4,4),add_scalebar=True)
```
<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/plt3.png" width="300"> <br>

(6) Core plotting functions included in `featureplot`.<br>
(6.1) Plot a feature with customized colormap with `plot1feature` and plot a region of interest with the argument `subset`.
```
cmap_1 = mcolors.LinearSegmentedColormap.from_list('name1',["black", "lime"], N=256)
fig=MALDIpy.featureplot.plot1feature(tissue_obj, cmap = cmap_1, max_num=41000, min_num=21000, figsize = (4.5,5))
fig=MALDIpy.featureplot.plot1feature_subset(tissue_obj,mz_use,cmap = cmap_1, max_num=41000, min_num=21000, figsize = (5,2.9),subset=[95,185,35,175])
```
<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/featureplot1.png" width="200">            <img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/featureplot2.png" width="350"> <br>

(6.2) Plot two features at the same time with customized colormap with `plot2features` and plot a region of interest with the argument `subset`.
```
cmap_1 = mcolors.LinearSegmentedColormap.from_list('name1',["black", "lime"], N=256)
cmap_2 = mcolors.LinearSegmentedColormap.from_list('name2',["black", "magenta"], N=256) 
fig=MALDIpy.featureplot.plot2features(tissue_obj, feats = [mz_use1,mz_use2], cmap=[cmap_1,cmap_2],
                          max_num_1=41000, min_num_1=21000, max_num_2=50000, min_num_2=25000)
fig=MALDIpy.featureplot.plot2features_subset(tissue_obj, feats = [mz_use1,mz_use2], cmap=[cmap_1,cmap_2],
                          max_num_1=41000, min_num_1=21000, max_num_2=50000, min_num_2=25000, subset=[95,185,35,175])
```
<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/featureplot3.png" width="200">            <img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/featureplot4.png" width="350"> <br>

(7) Project your cell cluster of interest onto a group of tissue sections with `projection.project_cluster_in_groups`.<br>
First, extract X/Y coordinates from the anndata and add them to adata.obs with `projection.add_coords`
```
adata = MALDIpy.projection.add_coords(adata)
```
Then, visualize the cluster of interest across multiple samples.
```
group=['sample1','sample2','sample3','sample4','sample5','sample6']
fig = MALDIpy.projection.project_cluster_in_groups(adata, cluster_id='1', cluster_obs_name='leiden', 
                                group_list=group, group_obs_name='sample_id', cmap='Reds')
```
<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/project_cluster_in_groups.png" width="1000"> <br>

# 4. Related links

(1) <a href="https://en.wikibooks.org/wiki/Metabolomics/Analytical_Methods/Mass_Spectrometry/MALDI-MS">MALDI-MS (Matrix-assisted laser desorption/ionization mass spectrometry)</a><br>
(2) <a href="https://metaspace2020.eu/">METASPACE - cloud platform for spatial metabolomics</a><br>
(3) <a href="https://scanpy.readthedocs.io/en/stable/index.html">Scanpy - Single-Cell Analysis in Python</a><br>
(4) <a href="https://portals.broadinstitute.org/harmony/">Harmony - integrating multiple high-dimensional datasets</a><br>

# 5. Citation

Ongoing project
<br><br>
<img src="https://github.com/TheHumphreysLab/MALDIpy/blob/main/readme_fig/AI_paint.jpg" width="500"> <br>
Adobe Firefly AI painting with command `a robot analyzing tissue section metabolites with imaging mass spectrometry`.
**************


Find us on Twitter: 
<br/>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false">@HumphreysLab</a>
<br/><br/>
