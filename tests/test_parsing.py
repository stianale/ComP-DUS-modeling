# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

from chai_lab.data.parsing.input_validation import (
    constituents_of_modified_fasta,
    identify_potential_entity_types,
)
from chai_lab.data.parsing.structure.entity_type import EntityType

from .example_inputs import example_dna, example_ligands, example_proteins, example_rna


def test_simple_protein_fasta():
    parts = constituents_of_modified_fasta("RKDES")
    assert parts is not None
    assert all(x == y for x, y in zip(parts, ["R", "K", "D", "E", "S"]))


def test_modified_protein_fasta():
    parts = constituents_of_modified_fasta("(KCJ)(SEP)(PPN)(B3S)(BAL)(PPN)KX(NH2)")
    assert parts is not None
    expected = ["KCJ", "SEP", "PPN", "B3S", "BAL", "PPN", "K", "X", "NH2"]
    assert all(x == y for x, y in zip(parts, expected))


def test_rna_fasta():
    seq = "ACUGACG"
    parts = constituents_of_modified_fasta(seq)
    assert parts is not None
    assert all(x == y for x, y in zip(parts, seq))


def test_dna_fasta():
    seq = "ACGACTAGCAT"
    parts = constituents_of_modified_fasta(seq)
    assert parts is not None
    assert all(x == y for x, y in zip(parts, seq))


def test_parsing():
    for ligand in example_ligands:
        assert EntityType.LIGAND in identify_potential_entity_types(ligand)

    for protein in example_proteins:
        assert EntityType.PROTEIN in identify_potential_entity_types(protein)

    for dna in example_dna:
        assert EntityType.DNA in identify_potential_entity_types(dna)

    for rna in example_rna:
        assert EntityType.RNA in identify_potential_entity_types(rna)