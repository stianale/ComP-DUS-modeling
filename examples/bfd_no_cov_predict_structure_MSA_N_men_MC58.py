from pathlib import Path

import numpy as np
import torch

from chai_lab.chai1 import run_inference
from chai_lab.data.parsing.msas.aligned_pqt import _merge_files_in_directory

# We use fasta-like format for inputs.
# - each entity encodes protein, ligand, RNA or DNA
# - each entity is labeled with unique name;
# - ligands are encoded with SMILES; modified residues encoded like AAA(SEP)AAA

# Example given below, just modify it


example_fasta = """
>protein|name=N_men_MC58
EKAKINAVRAALLENAHFMEKFYLQNGRFKQTSTKWPSLPIKEAEGFCIRLNGIARGALDSKFMLKAVAIDKDKNPFIIKMNENLVTFICKKSASSCSDGLDYFKGNDKDCKLLK
>dna|name=AT_DUS
ATGCCGTCTGAA
>dna|name=AT_SUD
TTCAGACGGCAT
""".strip()

fasta_path = Path("/tmp/example.fasta")
fasta_path.write_text(example_fasta)

# Define the directory containing your .a3m files
directory = '/media/stian/hgst6tb/chai-lab_wukevin/chai-lab/bfd_cov_50_MSA_N_men_MC58/final_no_cov'

# Call the function to merge and convert the files
_merge_files_in_directory(directory)

msa_path = Path("/media/stian/hgst6tb/chai-lab_wukevin/chai-lab/bfd_cov_50_MSA_N_men_MC58/final_no_cov")

output_dir = Path("/tmp/outputs")

candidates = run_inference(
    fasta_file=fasta_path,
    output_dir=output_dir,
    msa_directory=msa_path,
    # 'default' setup
    num_trunk_recycles=20,
    num_diffn_timesteps=400,
    seed=42,
    device=torch.device("cuda:0"),
    use_esm_embeddings=True,
)

cif_paths = candidates.cif_paths
scores = [rd.aggregate_score for rd in candidates.ranking_data]


# Load pTM, ipTM, pLDDTs and clash scores for sample 2
scores = np.load(output_dir.joinpath("scores.model_idx_2.npz"))