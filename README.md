# NCR_peptides_analysis
R code and data inputs for Higgins et al. 2023. Plant-derived, nodule-specific cysteine rich peptides inhibit growth and psyllid acquisition of ‘Candidatus Liberibacter asiaticus’, the citrus Huanglongbing bacterium. https://doi.org/10.1101/2023.06.18.545457 

# Download code and data

```{sh}
git clone https://github.com/higginss1/NCR_peptides_analysis.git
cd NCR_peptides_analysis
mkdir figures
```

You can now:
- open the `NCR_peptides_analysis.R` file in Rstudio
- set the working directory to `NCR_peptides_analysis/`
- Use run all to run code in the file serially
    - The script will:
        -  read in the necessary inputs from the cloned subdirectories
        -  output figures to `NCR_peptide_analysis/figures/`
        -  output any text files or other data to correct subdirectories

# Description of directory contents

| File or folder | Description |
| -------------- | ----------- |
| `NCR_peptides_analysis.R` | R code needed to reproduce the analysis |
| `Lcre_growth_assay_182_NCR_peptides/` | R code inputs and raw data for 96-well plate growth assays with *L. crescens* strain BT-1 |
| `NCR_Frac_Fac_128/` | R code inputs and raw data for plate growth assays with a fractional factorial design |
| `Detached_leaf_assay_data_inputs/` | R code and raw data for the NCR peptide analysis using a detached leaf assay |

## Lcre_growth_assay_182_NCR_peptides/

| File or folder | Description |
| -------------- | ----------- |
| `NCR_Peptides_*` | Directories containing raw 96-well plate growth information (.xls) for NCR peptide inhibition of *Liberibacter crescens* BT-1 and a text file containing all raw data modified for input into R |
| `Biomatik_182_NCR_peptide_data_long.txt` | NCR peptide sequence information in long format |
| `Biomatik_182_NCR_peptide_data_wide.txt` | NCR peptide sequence information in wide format |
| `Top_47_NCR_peptide_inhibitors.txt` | NCR peptide ID, percent growth rate inhibition, and peptide sequence of top 47 NCR peptides identified from 96-well plate assays |
| `all_peptides_GRAVY_score_edit.txt` | Detailed metadata for the NCR peptides investigated in the present study. Include details on the UniProt/Genbank accession ID, peptide sequence, the 20 amino acid sequence of the predicted antimicrobial peptide region, number of cysteine residues, hydrophobicity, charge (at pH 7), etc. |

## NCR_Frac_Fac_128/

| File or folder | Description |
| -------------- | ----------- |
| `NCR_Peptides_*` | Directories containing raw 96-well plate growth information (.xls) for top 10 NCR peptides utilized in a resolution V fractional factorial experimental design to test for synergistic effects of NCR peptides on inhibition of *Liberibacter crescens* BT-1 and a text file containing all raw data modified for input into R |
| `Biomatik_182_NCR_peptide_data_long.txt` | NCR peptide sequence information in long format |
| `Biomatik_182_NCR_peptide_data_wide.txt` | NCR peptide sequence information in wide format |
 `all_peptides_GRAVY_score_edit.txt` | Detailed metadata for the NCR peptides investigated in the present study. Include details on the UniProt/Genbank accession ID, peptide sequence, the 20 amino acid sequence of the predicted antimicrobial peptide region, number of cysteine residues, hydrophobicity, charge (at pH 7), etc. |

## Detached_leaf_assay_data_inputs/

| File or folder | Description |
| -------------- | ----------- |
| `.xls files` | Raw qPCR or RT-qPCR data for detached leaf assays performed using top-performing NCR peptides. |
| `SH_NCR_Peptide_DL_DLA_total_2021_09_28_qPCR_sample_sheet_final.txt` | A tab-delimited file containing information on the qPCR and RT-qPCR assays performed. |
| `detach_leaf_NCR_peptides_1_2.rds` | An Rdata formatted file containing qPCR and RT-qPCR data for two NCR peptides utilized in detached leaf assays |
| `detach_leaf_NCR_peptides_3-8.rds` | An Rdata formatted file containing qPCR and RT-qPCR data for an additional six NCR peptides utilized in detached leaf assays |




