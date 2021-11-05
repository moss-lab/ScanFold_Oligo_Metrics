# ScanFold_Oligo_Metrics
This ReadMe file contains information on inputs needed to run the script and the data contained in the output files from the ScanFold_Oligo_Metrics.py script.

This script was written in Python 3.8 to parse through ScanFold output data in tiled 18nt windows to determine metrics associated with potential ASOs.

The input files required to run this script are the FinalPartners.txt file, the input fasta file of the gene sequence, and the ScanFold.log file

To run the script you will type the following: python simple_oligo_metrics.py  FinalPartners.txt  Sequence.fasta  ScanFold.log  Outfile  IntegerWindowSize SequenceID
Here the input files names you enter will be the exact name of your file, the Outfile will be a name you decide, IntegerWindowSize will be 18, and SequenceID is the Mul_RS..... name

The output file of the scripts reports the average of average z-scores, average MFEs, and average EDs for the tiled user defined window, the total number of base pairs in the ScanFold predicted model for each tiled window, the number of possible base pairs per nucleotide for all nucleotides in the tiled user defined window according the ScanFold-Scan results, the reverse complement of the sequence in the tiled user defined window, the range of nucleotides being analyzed in the tiled window, and the sequence ID of the input.


The columns of data in the output file are as follows for each user defined window: 

1. avgZscore	
2. avgMFE	
3. avgED	
4. #ofBPs	
5. pos1BP/nt	
6. pos2BP/nt	
7. pos3BP/nt	
8. pos4BP/nt	
9. pos5BP/nt	
10. pos6BP/nt	
11. pos7BP/nt	
12. pos8BP/nt	
13. pos9BP/nt	
14. pos10BP/nt	
15. pos11BP/nt	
16. pos12BP/nt	
17. pos13BP/nt	
18. pos14BP/nt	
19. pos15BP/nt	
20. pos16BP/nt	
21. pos17BP/nt	
22. pos18BP/nt		
23. Reverse Complement
24. ntRange
25. TranscriptID


Below are the definition of each value type.

avgzscore: Average of the average z-scores across the tile user defined window.

avgMFE: Average of the average minimum free energies across the tiled user defined window.

avgED: Average of the average ensemble diversities across the tiled user defined window.

#ofBP: The number of total base pairs in the tiled user defined window according the to final ScanFold predicted model.

pos#BP/nt: The number of possible base pairs each nucleotide in the user defined window was predicted to make according to ScanFold-Scan.

Reverse Complement: The reverse complement of the sequence in the tiled user defined window. Included for ease of identifying siRNA seed sequences.

ntRange: The range of nucleotides being analyzed in each window

SequenceID: Transcript ID associated with each data window. This is included for ease of analysis when multiple files are concatenated into one larger file.
