# MultiCRISPR-EGA

## Introduction
MultiCRISPR-EGA is a tool for designing multiplexed CRISPR guide RNA (gRNA) arrays using an Elitist Genetic Algorithm (EGA). It supports efficient and simultaneous targeting of multiple genes, making it valuable for complex genome editing tasks. The tool offers both a user-friendly GUI and a CLI, and is compatible with various CRISPR-Cas systems

Datasets of varying sizes for Escherichia coli are available in the '**Gene_Info**' folder, and the generated multiplexed gRNA sequences can be found in the '**Result_xxx**' folders.


## Preparation
```

conda create --name multicrisprEGA python==3.10

conda activate multicrisprEGA

conda install --file requirements.txt
```

## Additional Software Installation
Install ViennaRNA and Ghostscript

1. Download the ViennaRNA from [ViennaRNA official site](https://www.tbi.univie.ac.at/RNA/) and Ghostscript from [Ghostscript download page](https://www.ghostscript.com/download/gsdnld.html).

2. Follow the installation instructions provided in the downloaded package.

3. Add their bin directory to your system's environment variables.

We also have provided a pre-downloaded ViennaRNA package in the '**Additional Software and Package**' folder

### Dataset Preparation
To extract genes from a specific species and generate CRISPR target sequences, follow these steps:

#### 1. Extract Gene Data from NCBI
    • Use the **Extract_from_NCBI.ipynb** notebook to retrieve gene information from the NCBI database for the desired species. 

    • Save the extracted gene data as a CSV file in the **Gene_Info** directory, following the format of existing files (e.g., **Gene_Info/Gene_Info4/GeneInfo.csv**).

    • To identify potential off-target sites, place the corresponding genomic data in the **Gene_Info/** folder.

#### 2. Generate CRISPR Target Sequences (Spacers)
Once you have the '**your_dataset_file.csv file**', use the following command to generate CRISPR target sequences (spacers):
```
cd CRISPR_library

python spacer_extract.py --geneInfoFile [Gene_Info/../your_dataset_file.csv] --PAM [PAM sequences] --pamPosition [before|after] --targetLength [target length] --outputFile [Gene_Info/../your_spacer_file.csv]
```

## Run the Application
This tool can be operated in two different modes depending on your preference:

### GUI Mode
To use the graphical user interface, simply run the following command:

```
python MultiCRIPSR_EGA_GUI.py
```
This will launch the graphical user interface, allowing you to input parameters and execute tasks through a user-friendly interface.

### CLI Mode
For command-line usage, execute the following command:

Step 1: Calculate Efficiency Scores

``` 
cd Efficiency

python CalEffScore.py --spacerFile  [../../Gene_Info/../your_spacer_file.csv] --outputFile [../../Gene_Info/../efficent_score.csv]
``` 

Step 2: Identify Off-Target Sites
``` 
cd Off_target

python Find_off_target_sites.py --spacerFile  [../../Gene_Info/../your_spacer_file.csv]  --fastaFile [../../Gene_Info/your_genomic_data_file.fasta ] --offinderInputFile [../../Gene_Info/../offinder_input.txt] --matchPattern [] --DNA_bumps_num [] --RNA_bumps_num [] --mismatches [] --offinderOutputFile [../../Gene_Info/../offinder_output.txt]
``` 

Step 3: Calculate Off-Target Scores
``` 
cd Off_target

python CalOffScore.py --offinderOutputFile  [../../Gene_Info/../offinder_output.txt] --outputFile [../../Gene_Info/../off_target_score.csv]
``` 

Step 4: Calculate Energy Scores
``` 
cd Energy

python CalFreeEnergy.py --spacerFile [../../Gene_Info/../your_spacer_file.csv] --outputFile [../../Gene_Info/../free_energy.csv]
``` 

Step 5: Combine All Scores into a CSV File
``` 
python spacerInfo_combine.py --spacerFile  [../Gene_Info/../your_spacer_file.csv] --off_target_score_file [../Gene_Info/../off_target_score.csv] --eff_score_file [../Gene_Info/../efficent_score.csv]  --energy_file [../Gene_Info/../free_energy.csv] --combineInfo_file_path  [../Gene_Info/../spacer_info_combine.csv]
``` 

Step 6: Filter Spacers Based on Score Ranking
``` 
python spacer_filter.py --ori_spacer_file [../Gene_Info/../your_spacer_file.csv] --combineInfo_file [../Gene_Info/../spacer_info_combine.csv] --filtered_file [../Gene_Info/../spacer_filtered.csv] --cleaned_spacer_file [../Gene_Info/../spacer_cleaned.csv] --offtarget_threshold 50 --efficiency_threshold 10 --energy_threshold 10
``` 

Final: Running EGA
``` 
python MultiCRIPSR_EGA_CLI.py --insulatorSeq 'AGCUGUCACCGGAUGTGCUUUCCGGUCUGAUGAGUCCGUGAGGACGAAACAGCCUCUACAAAUAAUUUUGUUUAA' --repeatSeq 'GUGUCAUAGCCCAGCUUGGCGGGCGAAGGCCAAGAC' --geneSpacerFile [Gene_Info/../spacer_cleaned.csv] --populationSize 100 --maxIterations 10 --crossover 0.7 --mutation 0.1 --tournamentSize 4
```

