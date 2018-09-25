# Size of each piece

Each piece is comprised of a number of pathway pairs. Their number should produce an input for improved_piece_PCxN_estimates01.R which won't stop the script due to memory issues. Also, since our aim is to run as many 01.Rs as we can concurrently, we shouldn't use the maximum memory for each piece.