# README of Correlated Protein Modules Revealing Functional Coordination of Interacting Proteins are Detected by Single-cell Proteomics

The proteomic datasets can be found at MassIVE (https://doi.org/doi:10.25345/C55717S21). 

The code for processing data can be found at GitHub (https://github.com/dionezhang/CPM).

## Description of the project 

We analyzed human K562 cells, 293T cells and mouse oocytes to measure the correlations between the translational levels of any pair of proteins in a single mammalian cell.

List of analysis files

- The raw files have been processed with Proteome Discoverer and subsequent analyses were performed in R. The code for processing data can be found at GitHub (https://github.com/dionezhang/CPM).
    - Correlation matrix of human K562 cells, 293T cells and mouse oocytes can be found in Supporting Dataset of the manuscript.
    - Code for processing data and generate figure are named with the related figure (e.g. Figure1.R for Figure 1).

## Description of data files 

The raw MS data and the results were deposited in MassIVE (ID: MSV000089625). The repository is structured as follows: 

- raw/
    - K562/
        - The folder contains 140 .raw files as acquired by the MS instrument. Each file contains LFQ data (hence 1 single cell per file, 70 cells in total). Files with “load” in their names are corresponding files in the sample loading process.
    - 293T/
        - The folder contains 192 .raw files as acquired by the MS instrument. Each file contains LFQ data (hence 1 single cell per file, 96 cells in total). Files with “load” in their names are corresponding files in the sample loading process.
    - oocyte/
        - The folder contains 308 .raw files as acquired by the MS instrument. Each file contains LFQ data (hence 1 single cell per file, 154 cells in total). Files with “load” in their names are corresponding files in the sample loading process.


