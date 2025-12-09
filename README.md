# The Cis-regulatory Code of Mouse Immune Cells

## Summary

Differential gene expression allows multicellular organisms to generate diverse cell types and tissues from the same genome. This process is governed by intricate gene-regulatory networks (GRNs) comprising thousands of transcriptional and post-transcriptional regulators. At the transcriptional level, cis-regulatory elements (CREs)—including promoters and enhancers—serve as binding sites for transcription factors (TFs), which recognize specific DNA motifs to modulate chromatin structure and transcriptional initiation. While promoters are located near transcription start sites (TSSs), enhancers act at distal loci, often looping to interact with promoters and fine-tune transcriptional output. The combinatorial activity of these elements drives the dynamic transcriptional programs essential for development and cell differentiation.

Integrating ATAC-seq with RNA-seq enables the mapping of chromatin accessibility to transcriptional output, yielding a mechanistic view of gene regulation [[1](#1),[2](#2)]. Such integrative datasets have illuminated regulatory programs in cancer, development, and immune differentiation (e.g., Corces et al. 2018; Cusanovich et al. 2018; [Yoshida et al. 2019](#3)). However, the high dimensionality and nonlinearity of these data demand advanced computational approaches. Deep learning models, as demonstrated by [*Maslova et al. (2020)*](#5), can capture complex dependencies between chromatin state, motif composition, and gene expression. By learning hierarchical representations directly from sequence and accessibility data, these models outperform traditional statistical frameworks in predicting regulatory activity and inferring GRN architecture, providing a powerful foundation for decoding transcriptional regulation at scale.

<div class="figure" style="text-align: center">
<img src="images/AI-TAC.png" alt="AI-TAC (from Maslova et al. 2020)" width="80%" />
<p class="caption"> <b>Fig. 2.</b> Using Deep Learning models to understand the cis-regulatory grammar of immune cells (Maslova et al. 2020) </p>
</div>



## Ideas for objectives and work plan

In this project, we aim to investigate the gene regulatory landscape of immune cells using paired ATAC-seq and RNA-seq measurements and analyze those with deep learning based sequence-to-function models. This study seeks to understand how CREs interplay to drive transcriptional changes during differentiation, and the motifs, their context, and their combinations that are involved.

To investigate this, you should try to look into one or several of the following approaches: 

1. Classical data analysis
2. Sequence-to-Activity model for ATAC-seq
   1. Train model and compare different modeling strategies
      1. From pre-processed data matrix or from more granular bigwigs
   2. Generate sequence attributions for pre-selected sets of sequences related to your scientific question
   3. Extract motifs, cluster, and determine transcription factors involved in the mechanism
   4. Confirm findings with RNA-seq
3. Sequence-to-Expression model RNA-seq
   1. Train model for RNA-seq
      1. Jointly train with ATAC-seq
      2. pre-train on ATAC-seq and use transfer learning
      3. Determine regulatory elements contributing gene expression changes
         1. Analyze their motifs, and identify transcription factors
      4. Derive cell lineage specific Gene Regulatory Networks
4. Compare Results from Sequence-to-Function models for Immune Cells between human and mice
   1. Do the same transcription factors regulate differential gene regulation?
   2. Build model for  [Calderon et al. 2019](#4) and compare motifs controlling differential gene expression


## Datasets

The data can be visualized in the UCSC genome browser, the link to these data can be found here: http://rstats.immgen.org/Chromatin/chromatin.html.

You can download the datasets for the analysis from here:

- [Processed ATAC-seq data and called peaks](https://sharehost.hms.harvard.edu/immgen/ImmGenATAC18_AllOCRsInfo.csv)
- [Processed RNA-seq data](https://www.cell.com/cms/10.1016/j.cell.2018.12.036/attachment/4392da81-c56e-471a-b1df-0e72f03ecd77/mmc2.csv)
- [Summary of Immune Cell Populations Profiled by ATAC-Seq and Their QC Matrices](https://www.cell.com/cms/10.1016/j.cell.2018.12.036/attachment/e5df7329-d77d-40b3-a03a-34bdbe4b402c/mmc1.xlsx)
- [Transcription Start Sites](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz)

	+ The columns in this file are 1) Gene Name , Transcript Name, Chromosome, Strand, 5' transcript Start, 3' Transcript Start, Coding Region Start, Coding Region End, Exon Count, Exon Starts, Exon Ends
	+ Use Strand information to determine TSS

- [Chromvar TF motif associations for all OCRs](https://sharehost.hms.harvard.edu/immgen/ImmGenATAC18_AllTFmotifsInOCRs.txt)
- [AI-TAC (Pre-trained model and software)](https://github.com/smaslova/AI-TAC)
- Processed [BigWig ATAC-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100738) files to train models genome-wide or visualize individual regions with IGV Browser. See [here](https://github.com/sasselab-teaching/GenomicS2F_seminar/blob/main/course_resources/FileDownloadGEO.md) for help with downloading them

## References

- **Hu et al., 2025**:  
  Hu, Y., Horlbeck, M. A., Zhang, R., et al. (2025). Multiscale footprints reveal the organization of cis-regulatory elements. *Nature, 638*, 779–786. [https://doi.org/10.1038/s41586-024-08443-4](https://doi.org/10.1038/s41586-024-08443-4) <a id="1">[1]</a>

- **Yan et al., 2020**:  
  Yan, F., Powell, D. R., Curtis, D. J., & Wong, N. C. (2020). From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis. *Genome Biology, 21*, 22. [https://doi.org/10.1186/s13059-020-1929-3](https://doi.org/10.1186/s13059-020-1929-3) <a id="2">[2]</a>

- **Yoshida et al., 2019**:  
  Yoshida, H., et al. (2019). The cis-Regulatory Atlas of the Mouse Immune System. *Cell, 176*(4), 897–912.e20. [https://doi.org/10.1016/j.cell.2018.12.036](https://doi.org/10.1016/j.cell.2018.12.036) <a id="3">[3]</a>

- **Calderon et al., 2019**:  
  Calderon, D., et al. (2019). Landscape of stimulation-responsive chromatin across diverse human immune cells. *Nature Genetics, 51*(10), 1494–1505. [https://doi.org/10.1038/s41588-019-0518-2](https://doi.org/10.1038/s41588-019-0518-2) <a id="4">[4]</a>

- **Maslova et al. 2020**:

  A. Maslova, R.N. Ramirez, K. Ma, H. Schmutz, C. Wang, C. Fox, B. Ng, C. Benoist, S. Mostafavi, & Immunological Genome Project, Deep learning of immune cell differentiation, Proc. Natl. Acad. Sci. U.S.A. 117 (41) 25655-25666, https://doi.org/10.1073/pnas.2011795117 (2020).<a id="5">[5]</a>

