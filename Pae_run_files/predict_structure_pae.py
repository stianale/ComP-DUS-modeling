from pathlib import Path
import numpy as np
import torch

from chai_lab.chai1 import run_inference

# Setup input and output paths
example_fasta = """
>protein|name=example-of-long-protein
AGSHSMRYFSTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQKYKRQAQTDRVSLRNLRGYYNQSEAGSHTLQWMFGCDLGPDGRLLRGYDQSAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGTCVEWLRRYLENGKETLQRAEHPKTHVTHHPVSDHEATLRCWALGFYPAEITLTWQWDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLTLRWEP
>protein|name=example-of-short-protein
AIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM
>protein|name=example-peptide
GAAL
>ligand|name=example-ligand-as-smiles
CCCCCCCCCCCCCC(=O)O
""".strip()

fasta_path = Path("/tmp/example.fasta")
fasta_path.write_text(example_fasta)

output_dir = Path("/tmp/outputs")

candidates = run_inference(
    fasta_file=fasta_path,
    output_dir=output_dir,
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    seed=42,
    device=torch.device("cuda:0"),
    use_esm_embeddings=True,
)

cif_paths = candidates.cif_paths
scores = [rd.aggregate_score for rd in candidates.ranking_data]


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
