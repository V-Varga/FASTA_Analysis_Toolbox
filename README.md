# fasta_an.py
Author: Virág Varga

Date: 28.10.2020

## Description

The program 'fasta_an.py' is built as a toolbox to perform various analyses on a DNA nucleotide FASTA file.

There are 3 types of analyses that version 1.0 of fasta_an.py can perform:

  - 6-Frame Translation [-t]: This function translates all sequences in the FASTA file into all 6 possible reading frames. It then evaluates which protein sequence(s) is/are the most likely to be the real reading frame(s), based on the lengths of the translated protein sequences.

  - Nucleotide Frequency Calculation [-c]: This function calculates the frequency of the 4 standard nucleotide bases present in each sequence in the FASTA file. The GC content of the file is also calculated. A pdf file containing pie charts displaying graphical representations of these frequencies is created.

  - BLAST @ NCBI [-b]: This function uses BioPython to BLAST the FASTA file against the NCBI database. Results are outputted to a user-specified file in XML format. It is therefore recommended that this file have a .xml extension.

In addition, the program has a quality-checking function integrated into the 3 functions mentioned above. This ensures that the input FASTA file is of the proper quality.

## Usage

Prior to use, 2 external python packages must be installed:
 - matplotlib (https://matplotlib.org/)
 - BioPython (https://biopython.org/)

Both the matplotlib library and the BioPython tool can be installed via `conda`. For further installation instructions and more information on their functionality and use, please visit their respective websites (linked above). 


This program can be run from the command line in the following ways:

```./fasta_an.py [-h] [-t] [-c] [-b] [-v] input_file output_file```

OR

```python fasta_an.py [-h] [-t] [-c] [-b] [-v] input_file output_file```

![program help page](https://github.com/V-Varga/FASTA_Analysis_Toolbox/blob/main/FASTA_an_imgs/fasta_an_help.png)

## Results/Visuals

A few images showing sample results are provided below:

Results of the 6-frame translation:
![results of fasta_an.py -t](https://github.com/V-Varga/FASTA_Analysis_Toolbox/blob/main/FASTA_an_imgs/results_t.png)

This portion of the program translates the sequences from the input DNA nucleotide FASTA into the protein products of the 6 possible reading frames. It also highlights which sequence(s) is/are most likely to be the real reading frame, based on which protein sequence(s) is are the longest.

Results of the nucleotide frequency calculations:
![results of fasta_an.py -c](https://github.com/V-Varga/FASTA_Analysis_Toolbox/blob/main/FASTA_an_imgs/results_c.png)

![results of fasta_an.py -c](https://github.com/V-Varga/FASTA_Analysis_Toolbox/blob/main/FASTA_an_imgs/results_c_graphics.png)

This portion of the program calculates the frequency with which each of the 2 standard DNA nucleotides appear in each sequence in the FASTA file, as well as the GC content of each sequence. A companion pdf will be made to the user-specified output file that includes graphical depictions of the relative nucleotide frequencies for each sequence.

Results of the BLAST:
![results of fasta_an.py -b](https://github.com/V-Varga/FASTA_Analysis_Toolbox/blob/main/FASTA_an_imgs/results_b.png)

This portion of the program uses BioPython to perform a BLASTn analysis of the input FASTA sequence, the results of which are outputted in XML format. It is therefore suggested that the specified output file have a '.xml' extension.

## Support

With questions, please contact the author of this program at:

virag.varga.bioinfo@gmail.com

## Project Status

Current Version: 1.0

While I have no explicit plans for a timeline to continue development of this program, I may eventually include other functions I find useful in the analysis of nucleotide FASTA files, or return to work out some of the known bugs and limitations.

## License

[MIT]
(https://choosealicense.com/licenses/mit/)
