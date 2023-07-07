import pynmrstar
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
import logging
import os
from multiprocessing import Pool, cpu_count
import numpy as np

logging.basicConfig(level=os.environ.get("LOGLEVEL", "WARNING"))


class Protein_Entry:
    def __init__(self, bmrb_id):
        self.bmrb_id = bmrb_id
        self.entry = None
        self.assignment = None
        self.orig_assignment = None
        self.anonymised_assignment = None
        self.sequence_from_bmrb = None
        self.sequence_from_assignment = None

        self._get_entry()

    def _get_entry(self):
        self.entry = pynmrstar.Entry.from_database(
            self.bmrb_id, convert_data_types=True
        )
        logging.info(f"Got Entry for {self.bmrb_id}")

    def get_assignment(
        self, tags=["Comp_index_ID", "Comp_ID", "Atom_ID", "Val", "Entry_ID"]
    ):
        data = self.entry.get_loops_by_category("Atom_chem_shift")[0]
        assignment_df = pd.DataFrame(data.get_tag(tags), columns=tags)

        # Convert 3 letter sequence to 1 letter sequence
        assignment_df["Comp_ID"] = assignment_df["Comp_ID"].map(str.capitalize)
        assignment_df["Comp_ID_sl"] = assignment_df["Comp_ID"].map(protein_letters_3to1)

        assignment_df = assignment_df.astype(
            {
                "Comp_index_ID": "int32",
                "Comp_ID": "str",
                "Atom_ID": "str",
                "Val": "float32",
                "Entry_ID": "int32",
                "Comp_ID_sl": "str",
            }
        )

        self.assignment = assignment_df
        self.orig_assignment = assignment_df.copy()
        logging.info(f"Got Assignment for {self.bmrb_id}")

    def get_sequence_from_bmrb(self):
        entity_saveframe = self.entry.get_saveframes_by_category("entity")[0]
        self.sequence_from_bmrb = "".join(
            entity_saveframe["Polymer_seq_one_letter_code"][0].split()
        )
        return self.sequence_from_bmrb

    def get_sequence_from_assignment(self):
        if self.assignment is None:
            raise TypeError(
                "Assignment is None. This method must be run after get_assignment(). "
            )
        amino_acid_list = self.assignment.groupby(["Comp_index_ID"])[
            "Comp_ID_sl"
        ].first()
        self.sequence_from_assignment = "".join(amino_acid_list.to_list())
        return self.sequence_from_assignment

    def remove_unwanted_atom_types_from_assignment(
        self, atom_list=["N", "H", "C", "CA", "CB"]
    ):
        self.assignment = self.assignment[self.assignment["Atom_ID"].isin(atom_list)]

    def get_anonymised_assignment(self):
        """
        Steps to anonymise Assignment

        1. Extract out N, H, Ca,Cb and C atoms
        2. For residue 2 onwards, extract the chemical shifts for N, H, Ca, Cb and C for that residue. Also extract Ca, Cb and C for the previous residue
        3. Rename Ca, Cb and C from the previous residue as pCa, pCb and Pc
        """
        df_list = []

        residues = sorted(self.assignment["Comp_index_ID"].unique())
        for i, residue in enumerate(residues):
            if i == 0:
                pass
            else:
                tmp_df = (self.assignment["Comp_index_ID"] == residue) | (
                    (self.assignment["Comp_index_ID"] == residues[i - 1])
                    & (self.assignment["Atom_ID"].isin(["CA", "CB", "C"]))
                )

                tmp_df1 = self.assignment[tmp_df]
                mask = (tmp_df1["Comp_index_ID"] == residues[i - 1]) & (
                    tmp_df1["Atom_ID"].isin(["CA", "CB", "C"])
                )
                tmp_df1.loc[mask, "Atom_ID"] = "p" + tmp_df1.loc[mask, "Atom_ID"]
                tmp_df1.loc[mask, "Comp_index_ID"] = residue

                df_list.append(tmp_df1)

        self.anonymised_assignment = pd.concat(df_list)

    def encode_atom_id(self):
        atom_id_encode_dict = {"H": 1, "N": 2, "C": 3, "CA": 4, "CB": 5}
        self.assignment["Atom_ID_encoded"] = self.assignment["Atom_ID"].map(
            atom_id_encode_dict
        )

    def get_dummy_comp_id(self):
        # create a dummy residue variable
        dummy_residue_values = np.arange(
            0, len(self.anonymised_assignment["Comp_index_ID"].unique())
        )
        np.random.shuffle(dummy_residue_values)

        mapping_dict = dict(
            zip(
                self.anonymised_assignment["Comp_index_ID"].unique(),
                dummy_residue_values,
            )
        )

        self.anonymised_assignment["Dummy_Comp_index_ID"] = self.anonymised_assignment[
            "Comp_index_ID"
        ].map(mapping_dict)


def process_protein(protein_id):
    print(f"{protein_id}")

    # Load Protein entry
    protein = Protein_Entry(protein_id)

    # Get assignment and remove unwanted atom types. Handle error if assignment doesn't exist
    try:
        protein.get_assignment()
        protein.remove_unwanted_atom_types_from_assignment()

        # Get Sequences and print out to dataframes
        protein.assignment["sequence_from_bmrb"] = protein.get_sequence_from_bmrb()
        protein.assignment[
            "sequence_from_assignment"
        ] = protein.get_sequence_from_assignment()

        # Calculate anonymised assignment and assign a dummy residue value
        protein.get_anonymised_assignment()
        protein.get_dummy_comp_id()

        return protein
    except IndexError as e:
        print(f'{e}: This likely means an assignment doesn"t exist. Skipping entry. ')
        return None


def build_assignments_csvs_parallel():
    df = pd.read_csv(
        "comparable_af2_nmr_structures.txt",
        delim_whitespace=True,
        names=[
            "UniProt-ID",
            "PDB-ID",
            "chain",
            "BMRB-ID",
            "residue-range-AF2",
            "residue-range-NMR",
        ],
        skiprows=1,
    )

    # Initialise Lists
    all_assignments_list = []
    all_anonymised_assignments_list = []

    n_cores = cpu_count() - 1

    with Pool(processes=n_cores) as pool:
        for protein in pool.imap(process_protein, df["BMRB-ID"]):
            if protein is None:
                pass
            else:
                all_assignments_list.append(protein.assignment)
                all_anonymised_assignments_list.append(protein.anonymised_assignment)

    # Create a single dataframe and write to csv for both assignment and anonymised assignment
    all_assignments_df = pd.concat(all_assignments_list)

    # Reorder column order and save to disk
    all_assignments_df = all_assignments_df[
        [
            "Entry_ID",
            "Comp_index_ID",
            "Comp_ID",
            "Comp_ID_sl",
            "Atom_ID",
            "Val",
            "sequence_from_assignment",
            "sequence_from_bmrb",
        ]
    ]
    all_assignments_df.to_csv("All_Assignments.csv")

    # Repeat for Anonymised Assignment
    all_anonymised_assignments_df = pd.concat(all_anonymised_assignments_list)

    # Reorder column order and save to disk
    all_anonymised_assignments_df = all_anonymised_assignments_df[
        [
            "Entry_ID",
            "Comp_index_ID",
            "Dummy_Comp_index_ID",
            "Comp_ID",
            "Comp_ID_sl",
            "Atom_ID",
            "Val",
            "sequence_from_assignment",
            "sequence_from_bmrb",
        ]
    ]
    all_anonymised_assignments_df.to_csv("All_Anonymised_Assignments.csv")


if __name__ == "__main__":
    build_assignments_csvs_parallel()
