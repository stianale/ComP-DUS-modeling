# Chai-1

Clone this repo, and then, in another directory, clone the https://github.com/chaidiscovery/chai-lab repo and install following the instructions in README.

Copy over the folders present here:

cp -r path/to/ComP-DUS-modeling/Chai path/to/chai-lab

cd path/to/chai-lab/Chai

The scripts to run ComP-DUS inferences are found in ComP_scripts. Use run_script.bash to run the inferences, ie. for the protein and DNA specified in the python config file paired_ComP_DUS_cov_50_paired_bfd_cov_50_MSA_B_den_ComP.py:

```bash run_script.bash ComP_script/paired_ComP_DUS_cov_50_paired_bfd_cov_50_MSA_B_den_ComP.py```

This will run 10 replicates of the inference specified in the python config file provided as the only argument to running this bash script.

# AF3

For modeling with AF3, go to https://alphafoldserver.com

# RF2NA and HADDOCK

In another directory, clone the repo https://github.com/haddocking/haddock3 and install following the instructions given.

In yet another directory, clone the GitHub repo https://github.com/uw-ipd/RoseTTAFold2NA and install following the instructions given.

```cd HADDOCK```

The RF2NA + HADDOCK pipeline is run using the script RF2NA_HADDOCK.bash. Pre-requisites are PyMOL (https://github.com/schrodinger/pymol-open-source), the programming-language R (https://www.r-project.org) with bio3d installed (https://cran.r-project.org/web/packages/bio3d/index.html), reduce (https://github.com/rlabduke/reduce) and pdb-tools (https://github.com/haddocking/pdb-tools).

You need to replace the RF2NA paths (/path/to) in the RF2NA_HADDOCK.bash script to your specific location. 


To run RF2NA + HADDOCK for B. denitrificans ComP (B_denitrificans_ComP.fa) modelled with AG-DUS (AG_DUS.fa), using the HADDOCK ambiguous restraints file B_denitrificans_ComP_restraints.list:

```bash RF2NA_HADDOCK.bash B_denitrificans_ComP AG_DUS B_denitrificans_ComP_restraints.list```

This will:
1. Run 100 RF2NA replicates.
2. Choose the best RF2NA output model (lowest PAE value)
3. Store the RF2NA-predicted ComP and DUS to separate PDB files
4. Run normal mode analysis on both the protein and DNA predicted by RF2NA
5. Make ensemble PDB files, containing normal mode conformations.
6. Use HADDOCK to cross-dock all sampled normal modes of ComP and DUS

You can use the korn shell script 200_best_models.ksh to store the 200 best HADDOCK models (in terms of HADDOCK score) to an ensemble PDB file:

ksh93 200_best_models.ksh B_denitrificans_ComP

# DockQ - internal consistency and cross-platform comparisons

Pairwise DockQ was calculated for all top models:

1. Within the same platform
2. Across different platforms,

using the flag --mapping ABC:ABC. HADDOCK models were not included in DockQ calculations.
