# Predicting the translation efficiency of messenger RNA in mammalian cells

## Summary
Protein abundance in mammalian cells is regulated at multiple levels, including transcription, mRNA stability, translation, and protein degradation. Among these layers, translation efficiency, which is a measure of how effectively an mRNA molecule is converted into protein, is one of the major determinants of protein abundance. The paper “Predicting the translation efficiency of messenger RNA in mammalian cells” (Zheng et al., 2025) explores how mRNA encoded features, such as sequence composition, codon usage, untranslated region (UTR) elements, and RNA structure, affect how efficiently an mRNA is translated into protein. Using computational models, they identify some of the sequence and structural determinants that drive translation efficiency, allowing for the prediction of protein output directly from mRNA features.

<div class="figure" style="text-align: center">
<img src="images/RiboNN.png" alt="RiboNN (from Zheng et al., 2025)" width="80%" />
<p class="caption"> <b>Fig. 2.</b>  Schematic representation showing the workflow of transcriptome-wide TE calculations for the human and mouse and TE correlations in humans. (Maslova et al. 2020) </p>
</div>

## Ideas for objectives and work plan

In this project, we aim to investigate the sequence determinants of translation regulation in human and mouse cell lines using matched Ribo-seq and RNA-seq measurements and analyze those with classical ML/deep learning based sequence-to-function models:

To investigate this, you should try to look into one or several of the following approaches: 

1. Use classical ML models (for instance: lasso, elastic net, random forest and LGBM regression model) to predict translation efficiency from UTR/transcript features: 5'UTR length, GC content, Kozak strength, uORFs. You can go beyond the features that were used in the paper and explore different algorithms. 

2. Use deep learning on raw mRNA sequence. Here, you can train a CNN or Transformer that takes only sequence input and predicts translation efficiency. In the paper, they stick to convolutional neural network, so exploring and comparing the perfromance of transformer-based models on similar tasks would be interesting.

3. Integrating RNA secondary structure into ML models. In the paper, the authors use RNAfold-derived structural features such as MFE to train classical models. It would be interesting to test whether adding information on paired and free bases as a sixth channel would improve the predictions. 

4. Transfer learning. The authors train the model on human and mouse data. It is also a possibility to perform transfer learning/fine-tune the model on data from other species. 

5. Motif discovery via neural networks. Use explainable AI methods to identify sequence motifs associated with high or low translation efficiency.

## Datasets

You can download the datasets for the analysis from here:

- Feature sizes, sequences, CV folds and TEs of human genes: https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-025-02712-x/MediaObjects/41587_2025_2712_MOESM3_ESM.xlsx

- Feature sizes, sequences, CV folds and TEs of mouse genes: https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-025-02712-x/MediaObjects/41587_2025_2712_MOESM4_ESM.xlsx

- Feature sizes, sequences, CV folds and TEs predicted by the human RiboNN models: https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-025-02712-x/MediaObjects/41587_2025_2712_MOESM5_ESM.xlsx

- Feature sizes, sequences, CV folds and TEs predicted by the mouse RiboNN models: https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-025-02712-x/MediaObjects/41587_2025_2712_MOESM6_ESM.xlsx


## References

Zheng, D., Persyn, L., Wang, J. et al. Predicting the translation efficiency of messenger RNA in mammalian cells. Nat Biotechnol (2025). https://doi.org/10.1038/s41587-025-02712-x

Agarwal V, Kelley DR. The genetic and biochemical determinants of mRNA degradation rates in mammals. Genome Biol. 2022 Nov 23;23(1):245. doi: 10.1186/s13059-022-02811-x. PMID: 36419176; PMCID: PMC9684954.

