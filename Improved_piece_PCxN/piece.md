# Improved_piece_PCxN
## Differences from the original Improved PCxN

1. Create matrix ID_map which contains the tissue numbers and GSEs that correspond to the task ID given by the job array

2. Changed input parameters. Replaced the tissue (first parameter) with id which corresponds to the task ID given by the job array

3. Using ID to extract tissue, GSE to be used in the script

4. Each improved_PCxN_estimates01_piece.R run takes care of a single GSE from a specific tissue 