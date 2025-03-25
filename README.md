# Gut_Phage
Gut Phage are fun





### Conda environment

Python 3.12 environment detailed in environment.yml file. Pypangraph installed via pip, it's useful for parsing Pangraph JSON files.



### File structure


`config.py` = specifies data, scripts, and analysis paths
`data_utils.py` = basic functions for parsing metadata into Python dictionaries, reformatting data, etc
`stats_utils.py` = functions for statistical analyses used in other scripts
`plot_utils.py` = functions for visualizing data. Useful to ensure same color schemes, heatmap styles, etc across plots.

```
.
├── scripts/
│   ├── plot_utils.py
│   ├── data_utils.py
│   ├── stats_utils.py
│   ├── config.py
│   ├── analyze_and_plot_intraspecies_structure.ipynb
│   └── etc..
├── data/
│   ├── Single_vOTUs/
│   │   ├── Core_Alignments/
│   │   │   ├── core_alignment_combined_sequences_vOTU-000003_vOTU-000139.fna
│   │   │   └── etc....
│   │   └── Pangraphs/
│   │       ├── pangraph_combined_sequences_vOTU-000001_vOTU-000075
│   │       └── etc....
│   ├── Pairs_vOTUs/
│   │   ├── Core_Alignments/
│   │   │   ├── core_alignment_vOTU-000001_n1020.fna
│   │   │   └── etc...
│   │   └── Pangraphs/
│   │       ├── pangraph_vOTU-000001_n1020.json
│   │       └── etc....
│   ├── votu_mean_lengths.txt
│   ├── votu_genome_counts.txt
│   └── etc....
├── analyses/
│   ├── fig1.png
│   ├── fig2.png
│   └── etc...
├── config.py
├── environment.yml
└── .gitignore
```
