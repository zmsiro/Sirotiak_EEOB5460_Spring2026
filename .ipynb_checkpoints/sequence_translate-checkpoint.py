######################## BCB 5460: Python Assignment ########################

# Your Mission: Complete and Document a Python Workflow in a Jupyter Notebook #

#####----- Overview ------#####

# Your colleague, Dr. X, has provided a partially written Python script to analyze biological data. As their bioinformatician, they have asked you to complete the work. 

# The analyses they need include:
# a. Translate cytochrome-b DNA sequences into amino acid sequences  
# b. Calculate GC content and amino acid molecular weight  
# c. Add these values to a pandas DataFrame containing penguin body mass  
# d. Create visualizations to explore relationships in the data  

# You will complete missing code, improve documentation, and produce a fully reproducible analysis.


#####----- Tasks -----#####
###--- Part 1: Functions ---###

# 1. Document Dr. X’s function  
#    - Explain what it does, its inputs, and outputs (in markdown and/or docstrings)

# 2. Write a translation function (using a loop) based on Dr. X's pseudo-code suggesti# on  
#    - Translate DNA → amino acids  
#    - Iterate over sequence in codons (groups of 3)  
#    - Use the Vertebrate Mitochondrial codon table  
#    - Ignore stop codons at the end  

# 3. Write an alternative translation function  
#    - Hint: BioPython has built-in tools for this  
 
# 4. Write a function that calculates the molecular weight of each amino acid sequence  
#    - Input: amino acid sequence  
#    - Output: molecular weight (float)  
#    - Consider `Bio.SeqUtils.ProtParam`
 
# 5. Compute GC content  
#    - Input: DNA sequence  
#    - Output: GC proportion (0–1)


###--- Part 2: Main Workflow ---###

# 6. Add two new columns to the penguin DataFrame: (1) molecular weight and (2) GC content.  
 
# 7. Process sequences and fill DataFrame  
#    - For each species:
#      - translate DNA → amino acids  
#      - compute molecular weight  
#      - compute GC content  
 
# 8. Plot body mass by species  
#    - Create a bar chart  
#    - Include labeled axes and a title  
#    - Answer in markdown:
#      - What is the smallest penguin species?  
#      - What is its geographic range?  
 
# 9. Plot molecular weight vs GC content  
#    - Create a scatter plot  
#    - x-axis: GC content  
#    - y-axis: molecular weight  
#    - Include labeled axes  
 
# 10. Export results  
#    - Save the final DataFrame as:  
#      `penguins_mass_cytb.csv`
 
# 11. Bonus (optional)  
#    - Add an additional analysis, visualization, or function  
#    - (+0.5 points if total score < 15)


###--- Repository & Submission ---###

# You must complete this assignment in a new GitHub repository.
 
# - Repository name:  
#   `LastName_EEOB5460_Spring2026`
 
# - Your repository must:
#   - contain only files for this assignment  
#   - include your Jupyter notebook and required data files  
#   - run from top to bottom without errors  
 
# - Submit the GitHub repository URL on Canvas by the deadline  


##-- Required Files --##

# Download and include:
# - `sequence_translate.py`  
# - `penguins_mass.csv`  
# - `penguins_cytb.fasta`  


#####----- Deliverables -----#####

# Your submission must include:
 
# - A Jupyter notebook that:
#   - runs without errors  
#   - includes clear markdown explanations  
#   - documents all functions and steps  
 
# - Completed code for all tasks (1–11)  
 
# - At least two plots:
#   - body mass by species  
#   - molecular weight vs GC content  
 
# - A final DataFrame containing:
#   - species  
#   - body mass  
#   - GC content  
#   - molecular weight  


##-- Grading Criteria --##
# - Code correctness (35%)  
# - Documentation and clarity (30%)  
# - Data analysis and visualization (20%)  
# - Reproducibility (10%)  
# - Following directions (5%)  
#   - Correct repository name  
#   - Clean, organized repo  
#   - Includes required files  
#   - Correct GitHub URL submitted  


##-- Disclaimer --##
# Not all required skills have been covered in class.  
# You are expected to use external resources (and cite them) to complete this assignment.


######################## Python Translate Script ########################

## Here's the start of our Python script. Thanks for completing it for me! - Dr. X
## IMPORTANT: install BioPython so that this will work

from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd

#%%%%%%%%%%%%%%%#
### FUNCTIONS ###
#%%%%%%%%%%%%%%%#

## 1 ##
## Dr. X: this gets sequences 
## Please properly document this function in the Jupyter notebook 
## Your descriptions of all functions should contain information about what the function does,
## as well as information about the return types and arguments.
def get_sequences_from_file(fasta_fn):
    sequence_data_dict = {}
    for record in SeqIO.parse(fasta_fn, "fasta"):
        description = record.description.split()
        species_name = description[1] + " " + description[2]
        sequence_data_dict[species_name] = record.seq
    return(sequence_data_dict)

## 2 ##
####### YOUR STRING-TRANSLATE FUNCTION ########
## Write a function that translates sequences
## All sequences start at codon position 1
## Complete a function that translates using a loop over the string of nucleotides
## Here is  some pseudo-code and suggestions
## feel free to change the function and variable names
# def translate_function(string_nucleotides): 
#     mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"] # this should work using BioPython (be sure to check what this returns)
#     for-loop through every 3rd position in string_nucleotides to get the codon using range subsets
#         # IMPORTANT: if the sequence has a stop codon at the end, you should leave it off
#         # this is how you can retrieve the amino acid: mito_table.forward_table[codon]
#         add the aa to aa_seq_string
#     return(aa_seq_string)

## 3 ##
####### YOUR ALTERNATIVE FUNCTION ########
## Is there a better way to write the translation function? (Hint: yes there is) 
## Perhaps using available BioPython library utilities?
## Please also write this function.


## 4 ##
####### YOUR COUNT AA ANALYSIS FUNCTION ########
## Write a function that calculates the molecular weight of each amino acid sequence.
## For this, you can use some BioPython functions. I think you can use the ProtParam module.
## For more info, check this out: http://biopython.org/wiki/ProtParam
## So you should import the following before defining your function:
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# def compute_molecular_weight(aa_seq):
#     # I think the ProtParam functions may require aa_seq to be a string.
#     # It may not work if the amino acid sequence has stop codons.
#     run the ProteinAnalysis() function on aa_seq
#	  return the molecular weight

## 5 ##
####### YOUR GC CONTENT ANALYSIS FUNCTION ########
## Write a function that calculates the GC-content (proportion of "G" and "C") of each DNA sequence and returns this value.


#%%%%%%%%%%%%%%#
###   MAIN   ###
#%%%%%%%%%%%%%%#

cytb_seqs = get_sequences_from_file("penguins_cytb.fasta") 

penguins_df = pd.read_csv("penguins_mass.csv") # Includes only data for body mass 
species_list = list(penguins_df.species)

## 6 ## 
## Add two new columns to the penguin DataFrame: (1) molecular weight and (2) GC content.
## Set the value to 'NaN' to indicate that these cells are currently empty.

## 7 ##
## Write a for-loop that translates each sequence and also gets molecular weight and computes the GC content
## of each translated sequence and adds those data to DataFrame
# for key, value in cytb_seqs.items():
#     aa_seq = nuc2aa_translate_function(value) # whichever function you prefer of #2 or #3
#     get the molecular weight of aa_seq
#     get the GC content of the DNA sequence
#     fill in empty cells in DF that you created above

## 8 ##
## Plot a bar-chart of the mass with the x-axes labeled with species names.
## *Q1* What is the smallest penguin species? 
## *Q2* What is the geographical range of this species?

## 9 ##
## Plot a visualization of the molecular weight (y-axis) as a function of GC-content (x-axis).

## 10 ##
## Save the new DataFrame to a file called "penguins_mass_cytb.csv"

## 11 - BONUS ##
## What else can we do with this dataset in Python? 
## Add functions or anything that might be interesting and fun. (optional)

