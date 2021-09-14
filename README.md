# mycoplasma_ribosome
## Analysis of _M. pneumoniae_ ribosome and polysome distributions inside cells and ribosomal protein extensions

Supporting scripts to a publication by 
Liang Xue, Swantje Lenz, Maria Zimmermann-Kogadeeva, Dimitry Tegunov, Patrick Cramer, Peer Bork, Juri Rappsilber, and Julia Mahamid

## Contents: 
### MatLab script _workflow_polysome_analysis.m_ 
This script analyses ribosome sequences in polysomes identified from _M. pneumoniae_ tomograms and performs
statistical analysis of ribosome pairs and comparison to randomly shuffled ribosomes. 
Each ribosome is classified based on its conformation according to the translation-elongation cycle. 

Required file in the **data_ribosomes** folder:

- motl_annoted_addpolysome_allcombined.txt

File contains coordinates of all identified ribosomes across _M. pneumoniae_ tomograms with class information. 
The format is described in the file _motl_format.ppt_. 

- polysome_tomoNum_pid_riboclass_sequence

File contains polysome sequences encoded with ribosome classes. The file format and classifications are explained in the file _info.ppt_. 

### python notebook _plot_ribosomal_protein_extensions_ (.ipynb and .html files)
This script makes a figure with schematic representation of amino acid extensions in _M. pneumoniae_ ribosomal proteins and their structural and functional information. 

Required file in the **data_rp_sequences** folder:

- rp_extensions20_allCN.fasta

Fasta file with amino acid sequences for ribosomal proteins with extensions > 20 amino acids as compared to _E. coli_ downloaded from NCBI (RefSeq)

- rp_extension_stats.csv

Extension table with information on extension lengths for the selected ribosomal proteins

- rp_xlink_information.csv

Information on cross links in ribosomal proteins from _O'Reilly & Xue et al, Science 2020_ (https://doi.org/10.1126/science.abb3758)

- rp_jpred_prediction.fasta

Fasta file with secondary structure predictions by JPRED (https://www.compbio.dundee.ac.uk/jpred/)

- rp_iupred_disorder_scores.csv

File with disorder scores for each amino acid in ribosomal proteins predicted by IUPRED (https://iupred2a.elte.hu/)

- rp_TNins_mapping.csv

Positions of transposon insertions in ribosomal proteins from _Miravet-Verde et al., Nucleic Acids Research 2020_ (https://doi.org/10.1093/nar/gkaa679)
