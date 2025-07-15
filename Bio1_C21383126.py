
"""
Program to identify and display all potential Open Reading Frames (ORFs) in a DNA sequence from a FASTA file

Outputs:
The frame number;,
The nucleotide and amino acid start and stop positions,
The length of the ORFs in nucleotide/ amino positions, 
The AA sequence for all potential Realistic O.R.F. 
"""

import sys
import tkinter
from tkinter import filedialog

orf_list = [] # Stores all ORF data

def translate_codons(nucleotides):
    """
    Translates a DNA sequence into an amino acid sequence
    """
    CodonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    AminoAcidSeq = "" 

    # Translate all groups of 3 nucleotides into amino acids
    for i in range(0, len(nucleotides), 3):
        if (i+2) < len(nucleotides): # Check if there is a group of 3
            codon = (f"{nucleotides[i]}{nucleotides[i+1]}{nucleotides[i+2]}")
            if codon in CodonTable:
                AminoAcid = CodonTable[codon] # translates the codon into an amino acid
                AminoAcidSeq += AminoAcid # Add to sequence

    return AminoAcidSeq
    
def find_orf(AminoAcidSeq, frame_num, nucleotides_len):
    """
    Finds ORFs from amino acid sequence of a specific frame
    """
    i = 0

    offsets = {1:1, 2:2, 3:3, -1:0, -2:1, -3:2} # Offset value for each frame
    offset = offsets[frame_num] # Get offset value for current frame

    # Go through each amino acid
    while i < len(AminoAcidSeq):
        foundStop = False
        if AminoAcidSeq[i] == "M": # Start codon
            start = i # Start of ORF
            for j in range(i, len(AminoAcidSeq)):
                if AminoAcidSeq[j] == "*": # Stop codon
                    foundStop = True
                    end = j # Don't include *
                    i = j # Continue on from after the stop
                    break
            else:
                # If no stop codon - go to the end
                end = len(AminoAcidSeq)

            if len(AminoAcidSeq[start:end]) > 10: # 10 Amino Acids = 30 nucleotides 
                if not foundStop: # No end codon
                    end -= 1 # move end back one - mark last aa as end
                aa_seq = AminoAcidSeq[start:end]
                aa_len = len(aa_seq) 

                # Get nucleotide positions
                if frame_num > 0: # Positive strand
                    nt_start = start * 3 + offset
                    nt_end = end * 3 + offset + 2
                else: # Reverse complement strand
                    nt_end = nucleotides_len - (end * 3 + offset + 2)
                    nt_start = nucleotides_len - (start * 3 + offset)

                nt_len = abs(nt_end - nt_start) + 1
                # Mark if no stop codon - same as NCBI online
                if not foundStop:
                    aa_end = f">{end - 1}"
                    nt_end = f">{nt_end}"
                else:
                    aa_end = end - 1

                # Save orf info
                orf_info = {
                    "frame": frame_num,
                    "aa_start": start,
                    "aa_end": aa_end,
                    "nt_start": nt_start,
                    "nt_end": nt_end,
                    "aa_len": aa_len,
                    "nt_len": nt_len,
                    "aa_seq": aa_seq
                }
                orf_list.append(orf_info)
        i += 1 # Move on to next amino acid

def reverse_complement(nucleotides):
    """
    Gets the reverse compliment of a DNA sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'} # Complement of each nucleotide
    reverse = ""
    for nucleotide in nucleotides:
        reverse += complement[nucleotide] # Get complement of each nucleotide
    return reverse[::-1] # Return the reverse complement


def main():

    root = tkinter.Tk()
    root.wm_withdraw() # This completely hides the root window

    # Use windows explorer to input the file name
    FileName = filedialog.askopenfilename(filetypes = [('All files','*.*')])
    root.destroy()

    # Open the file in reading text mode using error checking
    try:
        Fp1 = open(FileName,'r')

        Fp1.readline() # Skip descriptor line

        # Empty string for nucleotides
        nucleotides = ""

        # Remove end of line character from each line and add to one long string
        for line in Fp1.readlines():
            nucleotides += line.strip("\n")

        # Translate nucleotide codons into amino acids and find all possible ORFs

        # Forward frames
        find_orf(translate_codons(nucleotides), +1, len(nucleotides)) # Frame 1 - starting from 0
        find_orf(translate_codons(nucleotides[1:]), +2, len(nucleotides)) # Frame 2 - starting from 1
        find_orf(translate_codons(nucleotides[2:]), +3, len(nucleotides)) # Frame 3 - starting from 2

        # Reverse complement frames
        comp = reverse_complement(nucleotides)
        find_orf(translate_codons(comp), -1, len(nucleotides)) # Frame -1 
        find_orf(translate_codons(comp[1:]), -2, len(nucleotides)) # Frame -2
        find_orf(translate_codons(comp[2:]), -3, len(nucleotides)) # Frame -3

        print(f"{len(orf_list)} ORFs found") # How many found

        for idx, orf in enumerate(orf_list):
            print(f"ORF{idx+1}:")
            print(f"Frame: {orf['frame']}")
            print(f"AA Start: {orf['aa_start']}    AA End: {orf['aa_end']}")
            print(f"NT Start: {orf['nt_start']}   NT End: {orf['nt_end']}")
            print(f"AA Length: {orf['aa_len']}  NT Length: {orf['nt_len']}")
            print(f"AA Sequence: {orf['aa_seq']}\n")

    except IOError:
        # Error with file
        print("error unable to read file or file does not exist!!!")
        print("Exiting the program")
        stop = input()
        Fp1.close()
        sys.exit(1)

"""
TEST PLAN
    run the program and select the fasta file
    input the same file at https://www.ncbi.nlm.nih.gov/orffinder/
    check that the outputs match

"""     
#****************  executing the program ***************

main()
              