#!/bin/bash

# Define the lists of variables
DEMOG_LIST=("beis" "msmc" "constant")
REC_LABEL_LIST=("rec6" "rec8")
FREQ_LABEL_LIST=("denovo" "standingvariation1" "standingvariation5")
SEL_LABEL_LIST=("sel1" "sel2" "sel3" "sel4")

module load R

# Loop through all combinations of the variables
for DEMOG in "${DEMOG_LIST[@]}"; do
  for REC_LABEL in "${REC_LABEL_LIST[@]}"; do
    for FREQ_LABEL in "${FREQ_LABEL_LIST[@]}"; do
      for SEL_LABEL in "${SEL_LABEL_LIST[@]}"; do
        # Define the output directory
        DIR="../runs/${DEMOG}_${REC_LABEL}_${FREQ_LABEL}_${SEL_LABEL}"
        
        # Create the directory and navigate into it
        mkdir -p "$DIR"
        cd "$DIR" || exit
        
        # Load the R module and run the R script
        Rscript ../../scripts/createSeeds.R "$DEMOG" "$REC_LABEL" "$FREQ_LABEL" "$SEL_LABEL"
        
        # Go back to the base directory
        cd - || exit
      done
    done
  done
done
