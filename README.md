# MALDIpy: Single-cell analysis of MALDI-MS imaging mass spectrometry data
`MALDIpy` is a Python package for analyzing MALDI-MS imaging mass spectrometry data at the single-cell level.

Its function includes: metabolite species visualization, Scanpy-based single-cell analysis with UMAP clustering, projection of cell type annotations, integrative visualization and analysis.

# 1. Installation

To install `MALDIpy`, run it in a terminal

```
pip install MALDIpy
```

# 2. Key functions

(1) Create a MALDIpy object with `msi_data`.

```
maldi_obj = msi_data(raw_file, scale=10)
```

(2) Visualize any metabolite of interest in a MALDIpy object with `plt`.

```
maldi_obj.plt(mz=741.530654593237, figsize = (6,5), smooth=False, pos = 'lower left', remove_hs = True, cmap = "magma_r")
```
<img src="https://github.com/HaikuoLi/MALDIpy/blob/master/readme_fig/plt1.png" width="300"> <br>

(3) MALDIpy to AnnData conversion with `to_adata`.

```
adata = maldi_obj.to_adata(add_meta=True, csv_file=raw_file)
```

(4) Efficient single-cell quality control, dimension reduction and clustering with `single_cell`. This includes Harmony-based batch effect correction when processing multiple tissue sections.<br>
<img src="https://github.com/HaikuoLi/MALDIpy/blob/master/readme_fig/plt2.png" width="300"> <br>

(5) Project single-cell cluster annotation onto the tissue section with `projection`.
```
MALDIpy.projection.umap_projection(adata, file_name=raw_file,pltcmap=adata.uns['leiden_colors'],
                                   figtitle='Leiden Cluster Projection',figdpi=150, fig_size=(4,4),add_scalebar=True)
```
<img src="https://github.com/HaikuoLi/MALDIpy/blob/master/readme_fig/plt3.png" width="300"> <br>

# 3. Tutorials

(1) Analyze one single dataset with high efficiency -- all analysis can be done in 5-10 minutes on a labtop:
https://github.com/TheHumphreysLab/MALDIpy/blob/main/vignette/kidney_cortex.ipynb

(2) Integrative analysis of multiple MALDI-MS datasets:
https://github.com/TheHumphreysLab/MALDIpy/blob/main/vignette/kidney_integration.ipynb
