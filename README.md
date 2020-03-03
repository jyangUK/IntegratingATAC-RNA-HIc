
[![INSERT YOUR GRAPHIC HERE](https://personalpages.manchester.ac.uk/staff/jing.yang/Data/MRC_figure1.png)]()

# Simultaneous analysis of open chromatin, promoter interactions and gene expression in stimulated T cells implicates causal genes for rheumatoid arthritis
 scripts to accompany the paper [**"Simultaneous analysis of open chromatin, promoter interactions and gene expression in stimulated T cells implicates causal genes for rheumatoid arthritis"**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/manuscript.pdf) authored by <i>Jing Yang, Amanda McGovern,  Paul Martin, Kate Duffus, Xiangyu Ge, Peyman Zarrineh, Andrew P Morris, Antony Adamson, Peter Fraser, Magnus Rattray & Stephen Eyre</i>. Scripts are based on R and presented in Jupyter notebook. 

## Table of contents:
- [ATACseq](#ATACseq)
- [RNAseq](#RNAseq)
- [CHiC](#CHiC)
- [HiC](#HiC)
- [Linking_CHiC_ATACseq_RNAseq](#Linking_CHiC_ATACseq_RNAseq)

## ATACseq
- [**ATACseq_calculate_LRandBIC.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/ATACseq/ATACseq_calculate_LRandBIC.ipynb): calculate Loglikelihood ratio (LR) and bayesian information criteria (BIC) for ATAC-seq data.
- [**ATACseq_heatmap_Fig_3a.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/ATACseq/ATACseq_heatmap_Fig_3a.ipynb): heatmap of ATAC-seq data based on the clustering results from ATACseq_clustering_Fig_3b.ipynb.
- [**ATACseq_clustering_Fig_3b.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/ATACseq/ATACseq_clustering_Fig_3b.ipynb): hierarchical Gaussian Processing (GP) clustering of dynamic ATAC-seq time course data.
- [**ATACseq_supplementary_Fig3.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/ATACseq/ATACseq_supplementary_Fig3.ipynb): Supplementary Fig3 for ATAC-seq QC.

## RNAseq
- [**RNAseq_correlation_exon_intron_check_supplementary_Fig1a-c.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/RNAseq/RNAseq_correlation_exon_intron_check_supplementary_Fig1a-c.ipynb): illustrate the correlations between exons and introns counts data from RNA-seq samples. Figures are shown in Supplementary Fig. 1a-c. 
- [**RNAseq_supplementary_Fig1d.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/RNAseq/RNAseq_supplementary_Fig1d.ipynb): PCA plots for RNA-seq samples. Figure is shown in Supplementary Fig. 1d.
- [**RNAseq_quality_check_withdatafromYe_supplementary_Fig1e-i.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/RNAseq/RNAseq_quality_check_withdatafromYe_supplementary_Fig1e-i.ipynb): compare the gene expression data from this study with the data from **"Ye, C. J., et al. "Intersection of population variation and autoimmunity genetics in human T cell activation." Science 345.6202 (2014): 1254665"**.

## CHiC
- [**CHiC_qualitycheck_supplementaryFig5a.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/CHiC/CHiC_qualitycheck_supplementaryFig5a.ipynb): compare CHi-C interactions with similar data from **"Burren, O.S. et al, Chromosome contacts in activated T cells identify autoimmune disease candidate genes, Genome Biology, 18 (165) (2017)"** . Venn diagram shown in Supplementary Fig. 5a
- [**CHiC_clustering_supplementary_FIg5b.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/CHiC/CHiC_clustering_supplementary_FIg5b.ipynb): clustering of CHi-C data and generate Supplementary Fig. 5b

## HiC
- [**HiC_interaction_matrices.ipynb**](https://github.com/jyanguk/IntegratingATAC-RNA-HiC/blob/master/HiC/HiC_interaction_matrices.ipynb): generate the interaction matrices for Hi-C data at different time. Upper triangular part of Hi-C data of chr1 is shown in Fig. 2a
- [**Fig_2b_2c_HiC_ABcompartment_correlations_plot.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/HiC/Fig_2b_2c_HiC_ABcompartment_correlations_plot.ipynb): compute the Stratum adjusted Correlation Coefficient (SCC) between Hi-C datasets and generate Fig. 2b. Compute the correlation coeffient between AB compartment data and generate Fig. 2c.
- [**TADs_percentage_plot_Supplementary_Fig.6.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/HiC/TADs_percentage_plot_Supplementary_Fig.6.ipynb): illustrate the percentage of the intersected TADs/A/B compartments between different samples.   
- [**ABcompartment_SNPs_overlap.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/HiC/ABcompartment_SNPs_overlap.ipynb): illustrate the overlap of rheumatoid arthritis SNPs and A/B compartments at different times.      
## Linking_CHiC_ATACseq_RNAseq
Directory for producing correlation density plots between linked CHi-C, ATAC-seq and RNA-seq data.
- [**plot_CHiC_ATACseq_RNAseq_connections_Fig4a.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/Linking_CHiC_ATACseq_RNAseq/plot_CHiC_ATACseq_RNAseq_connections_Fig4a.ipynb): generate the correlation density maps of linked CHi-C, ATAC-seq and RNA-seq data under varied distance ranges around promoters. p values from Wilcoxon tests are also shown in the plots. 
- [**plot_foldchange_Fig4bc.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/Linking_CHiC_ATACseq_RNAseq/plot_foldchange_Fig4bc-fromoriginallinked.ipynb): display the log2 fold change of ATAC-seq vs gene and CHi-C vs gene, as shown in Fig. 4b-c
- [**plot_supplementary_Fig_10.ipynb**](https://github.com/jyangUK/IntegratingATAC-RNA-HiC/blob/master/Linking_CHiC_ATACseq_RNAseq/plot_supplementary_Fig_10.ipynb): display dynamic patterns of CHi-C and ATAC-seq associated with different gene clusters. 

For more clarification, please feel free to contact Jing Yang via : Jing.Yang@manchester.ac.uk
