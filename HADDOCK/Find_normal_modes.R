#!/usr/bin/Rscript

library(bio3d)

args <- commandArgs(trailingOnly = TRUE)
inputfile <- args[1]

setwd("/media/stian/hgst6tb/OneDrive/DUS/PhD/All_Neis/Representative_genomes/Pipeline_HADDOCK/")

## Fetch structure
pdb <- read.pdb(inputfile)

## Calculate all-atom normal modes
modes.aa <- aanma(pdb, outmodes = 'noh')

## Calculate all-atom normal modes with RTB approximation
modes.aa.rtb <- aanma(pdb, outmodes='noh', rtb=TRUE)

## Print modes
print(modes.aa)
## Plot modes
plot(modes.aa)

## Extract the base name of the input file
input_base_name <- gsub("\\.pdb$", "", basename(inputfile))

## Visualize modes
for (mode in 7:16) {
  outfile <- file.path(getwd(), paste0("Normal_mode_", input_base_name, "_", mode, ".pdb"))
  mktrj(modes.aa, mode = mode, pdb = pdb, file = outfile)
}
