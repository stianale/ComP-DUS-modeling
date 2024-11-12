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
>protein|name=N_subflava_NJ9703_ComP_trimmed
RSANLRAAHAALLENARFMEQFYAKKGSFKLTSTKWPELPVKEAGGFCIRMSGQAKGILEGKFTLKAVALDREAEPRVLRLNESLTAVVCGKMKGKGSCTDGEEIFRGNDAECRPFTG
>dna|name=AG_DUS
AGGCCGTCTGAA
>dna|name=AG_SUD
TTCAGACGGCCT
""".strip()

fasta_path = Path("/tmp/example.fasta")
fasta_path.write_text(example_fasta)

msa_path = Path("/media/stian/hgst6tb/chai-lab/uniref_cov_50_MSA_N_subflava_NJ9703_ComP_trimmed/paired_ComP_DUS")

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

# Retrieve individual arrays and convert to NumPy
pae_scores = candidates.pae.numpy()
pde_scores = candidates.pde.numpy()
plddt_scores = candidates.plddt.numpy()

# Compute the mean for each array
pae_mean = np.mean(pae_scores)
pde_mean = np.mean(pde_scores)
plddt_mean = np.mean(plddt_scores)

# Save PAE, PDE, and pLDDT scores along with their means to a single .npz file
np.savez(output_dir / "scores_with_means.npz",
         pae=pae_scores, pde=pde_scores, plddt=plddt_scores,
         pae_mean=pae_mean, pde_mean=pde_mean, plddt_mean=plddt_mean)
