# STIMS.jl by Rohan Maddamsetti

## STIMS is a Simple Test to Infer Mode of Selection (STIMS), designed for metagenomic time series of evolution in hypermutator bacterial populations.

For details, see our preprint: Discovery of positive and purifying selection in metagenomic time series of hypermutator microbial populations, by Rohan Maddamsetti and Nkrumah A. Grant.

### Examples of command-line usage:

julia STIMS.jl SLiM-5000gen-OnePercent-Hypermutator.csv SLiM_geneIDs.csv SLiM_positive_module.csv -o SLiM-5000gen-OnePercent-Hypermutator-positive.pdf

julia STIMS.jl SLiM-5000gen-OnePercent-Hypermutator.csv SLiM_geneIDs.csv SLiM_neutral_module.csv -o SLiM-5000gen-OnePercent-Hypermutator-neutral.pdf

julia STIMS.jl SLiM-5000gen-OnePercent-Hypermutator.csv SLiM_geneIDs.csv SLiM_purifying_module.csv -o SLiM-5000gen-OnePercent-Hypermutator-purifying.pdf

