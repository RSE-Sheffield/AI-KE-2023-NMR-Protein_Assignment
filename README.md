# AI-KE-2023-NMR-Protein_Assignment

This repository contains a python script for preparing data for testing automated NMR protein assignment. The original scope of the 2 week project was to design a machine learning model that could automate the NMR protein assignment problem.

The repository contains 2 main files:

- *nmr_prot_assignment.py* - Python script that generates two csv files:
    - All_Assignments.csv
    - All_Anonymised_Assignments.csv
- *comparable_af2_nmr_structures.txt* - Text file containing NMR protein IDs that we want to extract from the BMRB database.

There is also an image of the chemical shift distributions for each atom across all proteins.
- *Atom_Chemical_Shift_distributions.png* -

## Run Python File

The requirements.txt file has all the python modules need to run the python file.

```
pip install requirements.txt
```

Then simply run the python file:

```
python nmr_prot_assignment.py
```

These are the steps nmr_prot_assignment.py performs:

1. Loops through each proteins and pulls data from [BMRB](https://bmrb.io/) using the IDs in comparable_af2_nmr_structures.txt
2. Each protein has the assignments extracted and remove unwanted atoms so only C, Ca, Cb, N and H remain.
3. Sequences are also extracted into the dataframe.
4. All Assignments are then joined together for each protein to make an All_Assignments.csv
5. This process is repeated but for each protein, an anonymised assignment is also created. The assignment is anonymised by copying C, Ca and Cb from the previous residue to current residue and changing the name to pC, pCa and pCb. A dummy comp index value is also produced.

## The built assignments.csv

The columns for the output csv files are:

- *Index*
- *Entry_ID* - BMRB entry id
- *Comp_index_ID* - Residue Number
- *Dummy_Comp_index* - Dummy Residue Number
- *Comp_ID* - Three letter Amino Acid
- *Comp_ID_sl* - Single Letter Amino Acid
- *Atom_ID* - Atom Name
- *Val* - Chemical Shift
- *sequence_from_assignment* - Sequence derived from assignment
- *sequence_from_bmrb* - Sequence from BMRB
