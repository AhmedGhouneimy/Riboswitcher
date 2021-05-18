#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
import os
import csv
from collections import OrderedDict
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from input_variables import *
####### Define some functions #########
def get_trigger_region(pos, seq, len_trigger, len_sensed):
  """
  returns the 200 trigger region and the sp and ep of the trigger binding sites 
  """
    #critical_pt = (len_trigger - len_sensed)//2
  critical_pt = len_trigger //2 #half point in the trigger
  trigger = ""
  #trigger_start_index
  if pos >= critical_pt and pos <= len(seq)-critical_pt: #if it's far away from both ends --> trigger seq is in the middle
      start = pos-critical_pt#- int(0.5*len_sensed)
      #new_pos = start
      end = pos+critical_pt#+ int(0.5*len_sensed)
      trigger = seq[start :  end]
      sp = critical_pt#-int(0.5*len_sensed) #sp start pos in the short trigger seq
      ep =critical_pt+len_sensed#+int(0.5*len_sensed) #ep end pos in the short trigger seq
      print ('condition1')
  elif pos <= critical_pt:
      #new_pos = pos
      start = 0
      end = len_trigger
      trigger = seq[start:end]
      sp = pos
      ep = pos + len_sensed
      trigger_start_index = pos
      print ('condition2')
  elif pos >= len(seq) - critical_pt:
      #new_pos = pos
      start = len(seq)-len_trigger
      end = len(seq)
      trigger = seq[start:]
      sp = critical_pt + (pos % critical_pt)#len_trigger-(pos % critical_pt)
      ep = sp + len_sensed
      trigger_start_index = pos
      print ('condition3')
  return trigger, sp, ep
 
def Hamming_Distance(str1, str2):
  print(len(str1) == len(str2))
  #print(len(str2))
  mismatches = 0
  for i in range (min(len(str1),len(str2))):
    if str1[i] != str2[i]:
      mismatches +=1
  print ('mismatches = ', str(mismatches))
  #print(min(len(str1),len(str2))-mismatches)/min(len(str1),len(str2))
  return (min(len(str1),len(str2))-mismatches)/min(len(str1),len(str2))

def nupack_energy(nupack_input):
    command = 'energy -material rna1995 -T 37 -multi ' + nupack_input
    #command = 'defect -material rna1995 -T 37 -multi -mfe '+ dot_structure + nupack_input
    print(command)
    a = os.system(command)
    return(a)
def nupack_pfunc(nupack_input):
    command = 'pfunc -material rna1995 -T 37 -multi ' + nupack_input
    #command = 'defect -material rna1995 -T 37 -multi -mfe '+ dot_structure + nupack_input
    print(command)
    a = os.system(command)
    return(a)
def nupack_defect(nupack_input, dot_structure):
    command = 'defect -T 37 -multi -material rna1995 -mfe ' + nupack_input
    #command = 'defect -material rna1995 -T 37 -multi -mfe '+ dot_structure + nupack_input
    print(command)
    a = os.system(command)
    return(a)

def nupack_pairs(nupack_input):
    '''This function receives a NUPACK input file and calculates the structure with the minimum free energy.
    It is important to remember that NUPACK input files must be named with the ".in" extension but this extension
    should not be included when calling NUPACK. For example: when executing NUPACK with a file called "example.in",
    we only need to say "example".
    '''
    command = 'pairs -material rna1995 -T 37 -multi ' + nupack_input
    print(command)
    a = os.system(command)
    return(a)
def parse_pairs(nupack_output,start_pos, end_pos):
    '''
    Takes start_pos: start of the triger binding site, end_pos: end of the trigger binding site
    This function parses the output by NUPACK and returns a dictionary with the summary of the results,
    that is, the sequence, the free energy, the secondary structure, and the list of paired bases.
    '''
    handle = open(nupack_output, 'r')
    out_dict = OrderedDict() ################
    
    
    # This counter will help me load the data
    counter = 0
    
    pii = [] # (base : p) probability that base i is single stranded
    unpaired_probability = []
    for line in handle:
        # Retrieve the sequence
        if line.startswith("%"):
          counter = 0
        # Look for the separator that will start giving us the numbers we want
        else:
          counter +=1 
        
        if counter == 2:
            # Then I know this is the length of the sequence.
            counter = counter + 1
            #print("line is: " + line)
            Nplusone = float(line) + 1 
   
        elif counter >= 3:
            # Then I am looking at the table of paired bases

            pii.append(tuple(line.strip().split('\t')))
    handle.close()

    sumpii = 0
    countsumpii = 0
    for base_probability in pii:#we can add a condition to only up to 15 bp probability
      if float(base_probability[1]) == Nplusone and float(base_probability[0]) >=start_pos and float(base_probability[0]) <= end_pos : #(N+1) then i know this is the probability that base is unpaired
        unpaired_probability.append(base_probability)
      #base_probability[0]: base number, base_probability[1] base paired with, base_probability[2] base pairing probability     
        sumpii += float(base_probability[2])
        countsumpii +=1 # to check if the len sensed  == the number of probabilities returned
    #print(unpaired_probability)
    #l_dna = sumpii / (end_pos-start_pos) #calculates the local single strandedness
    l_dna = sumpii / countsumpii
    #print("countsumpii" + str(countsumpii))
    #print( "end_pos-start_pos" + str(end_pos-start_pos))
    return (l_dna)
    # Return the dictionary
    #out_dict['unpaired_probability'] = unpaired_probability
    
    #return(unpaired_probability) #we can also return sumpii: the sum of all probab of all unpaired bases

def get_switch_recognition_seq(trigger, sequence_type, length_unpaired):
    ''' This function receives a target trigger sequence and the type of molecule
    and obtains the RNA trigger for it.
    '''
    # if sequence_type == 'RNA':
    #     trigger_seq = Seq(trigger, generic_rna)
    #     return(trigger_seq.back_transcribe().reverse_complement().transcribe())
    # elif sequence_type == 'DNA':
    #     trigger_seq = Seq(trigger, generic_dna)
    #     return(trigger_seq.reverse_complement().transcribe())
    
    trigger_seq = Seq(trigger, generic_rna)
    return(trigger_seq.back_transcribe().reverse_complement().transcribe())
"""    
def get_full_switch(recognition_sequence, reporter, length_unpaired):
    ''' This sequence receives the recognition sequence and the reporter gene to design the 
    full toehold switch. In the meantime, I will work with all my sequences in RNA.
    '''
    # Define the rest of the parts of the sequences
    part1 = 'GGG'
    part2 = recognition_sequence
    
    # Get the reverse complement of the trigger to close the loop.
    # Remember to add the start codon
    rev_comp = str(part2[length_unpaired:].reverse_complement())
    
    # Parts 3 and 5 should be complementary
    # Adjust their length so that the hairpin has 19 total paired bases with 3 weak pairs at the top (AUA - UAU pairs)
    # Get the number of needed bases
    needed = 16 - (len(recognition_sequence) - length_unpaired)
    
    part3 = 'AUA' 
    part4 = 'CAGAAACAGAGGAGA'
    part5 = 'UAU' 
    
    if needed <= 3:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'A'
        part6 = needed*'U'
        
        # Check how many bases of the complement of part2 should go before the AUG
        comp_pos = 3 - needed
        part6 = part6 + rev_comp[0:comp_pos] + 'AUG' + rev_comp[comp_pos+3:]
    elif needed <= 6:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'A'
        part6 = needed*'U'
        
        # Check where the complement of part2 should start
        comp_pos = 6 - needed
        part6 = part6 + 'AUG' + rev_comp[comp_pos:]
    elif needed > 6:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'A'
        part6 = 'UUU' + 'AUG' + (needed-6)*'U' + rev_comp
    
    # part7 = 'ACCUGGCGGCAGCGCAAAAG'
    part7_comp = Seq(reporter[1:5]).back_transcribe().reverse_complement().transcribe()
    part7 = 'ACCUGGCGGCAG' + str(part7_comp) + 'AAAG'
    
    # This is an example that fails because it generates a stop codon
    # part7 = 'ACUAAGCGGCAG' + str(part7_comp) + 'AAAG'
    
    # Test for stop codons
    stop_codons = ['UGA', 'UAA', 'UAG']
    test_region = part6[3:] + part7
    
    # Split the test region into codons to check one at a time
    test_region_codons = [test_region[i:i+3] for i in range(0, len(test_region), 3)]
                          
    # Loop through the list and check if there are any stop codons in this reading frame
    for codon in test_region_codons:
        # assert not codon in stop_codons, 'The generated switch contains a stop codon' 
        if codon in stop_codons:
            return('Stop')
                          
    part8 = reporter
    
    full_toehold = part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8
    
    return(full_toehold)
"""
def get_full_switch(recognition_sequence, reporter, length_unpaired):
    ''' This sequence receives the recognition sequence and the reporter gene to design the 
    full toehold switch. In the meantime, I will work with all my sequences in RNA.
    '''
    # Define the rest of the parts of the sequences
    #part1 = 'GGG'
    part1 = Upstream_UTR
    part2 = recognition_sequence

    
    # Get the reverse complement of the trigger to close the loop.
    # Remember to add the start codon
    rev_comp = str(part2[length_unpaired:].reverse_complement())
    # Parts 3 and 5 should be complementary
    # Adjust their length so that the hairpin has 19 total paired bases with 3 weak pairs at the top (AUA - UAU pairs)
    # Get the number of needed bases
    #Well this means that the step ubove the middle bulge is rationally designed
    needed = 16 - (len(recognition_sequence) - length_unpaired) #Leave as is bs check which class of toeholds is that 
    
    #part3 = 'AUA' #before the loop
    #part4 = "CUUUAACCGCCACCAUGUAAAG"#awh 'CAGAAACAGAGGAGA' #checked with nupack --> Linear structure
    #part4 = "CCUAAUACCGCCACCAUGAUUAGG"
    part4 = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    #part5 = 'UAU' #after the loop
    part6 = rev_comp #awh
    if len(part6)%3 != 0:
      print('there will be a frame shift because the reverse complement region (length_paired in input_variables.py) of the swich is not devisible by 3')
    print('recognition sequence: \n' + recognition_sequence + '\n len = ' + str(len(recognition_sequence)))
    print ('rev_comp: \n' + rev_comp  + '\n len = ' + str(len(rev_comp)))
    part7 = "AACCUGGCGGCAGCGCAAAAG" # awh Green et al's Linker (mRNA)
    #part7 = "CCUAAUACCGCCACCAUGAUUAGG"
    ''' #awh
    if needed <= 3:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'A'
        part6 = needed*'U' 
        
        # Check how many bases of the complement of part2 should go before the AUG
        comp_pos = 3 - needed
        part6 = part6 + rev_comp[0:comp_pos] + 'AUG' + rev_comp[comp_pos+3:]
    elif needed <= 6:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'A'
        part6 = needed*'U'
        
        # Check where the complement of part2 should start
        comp_pos = 6 - needed
        part6 = part6 + 'AUG' + rev_comp[comp_pos:]
    elif needed > 6:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'A'
        part6 = 'UUU' + 'AUG' + (needed-6)*'U' + rev_comp
    '''
    # part7 = 'ACCUGGCGGCAGCGCAAAAG'
    #awh part7_comp = Seq(reporter[1:5]).back_transcribe().reverse_complement().transcribe()
    #awh part7 = 'ACCUGGCGGCAG' + str(part7_comp) + 'AAAG'
    
    # This is an example that fails because it generates a stop codon
    # part7 = 'ACUAAGCGGCAG' + str(part7_comp) + 'AAAG'
    
    # Test for stop codons
    stop_codons = ['UGA', 'UAA', 'UAG']
    test_region = part6 + part7 #awh part6[3:] + part7
    
    # Split the test region into codons to check one at a time
    test_region_codons = [test_region[i:i+3] for i in range(0, len(test_region), 3)]
                          
    # Loop through the list and check if there are any stop codons in this reading frame
    for codon in test_region_codons:
        # assert not codon in stop_codons, 'The generated switch contains a stop codon' 
        if codon in stop_codons:
            return('Stop')
                          
    part8 = reporter
    
    #full_toehold = part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8
    full_toehold = part1 + part2 + part4 + part6 + part7 + part8

    return(full_toehold)
def prepare_nupack_input(sequence, output_file):
    '''This function writes file with the sequence and the specifications to use only one instance of it.
    '''
    handle = open(output_file, 'w')
    handle.write('1\n')
    handle.write(str(sequence) + '\n')
    handle.write('1')
    handle.close()

# Write a function to prepare the NUPACK input to test binding to the mRNA
def prepare_nupack_input_two(toehold, mRNA, output_file):
    '''This function prepares an input for NUPACK to test how the proposed toehold binds to the mRNA.
    '''
    handle = open(output_file, 'w')
    handle.write('2\n')
    handle.write(str(toehold) + '\n')
    handle.write(str(mRNA) + '\n')
    handle.write('1 2')
    handle.close()



def check_matches(start_pos, end_pos, pair_list):
    '''This function will generate the list of pairs of positions that should be paired and compare
    it to the list of matches predicted by NUPACK when simulating binding.
    '''
    list_targets = []
    # The starting base of the matching region
    check_match = 4
    for target in range(end_pos, start_pos - 1, -1):
        list_targets.append((str(check_match), str(target)))
        # Move on to the next check
        check_match = check_match + 1

    right_matches = 0
    total_matches = 0
    for pair in list_targets:
        total_matches = total_matches + 1
        if pair in pair_list:
            right_matches = right_matches + 1
    
    return(round(right_matches * 100 / total_matches, 2))


def nupack_mfe(nupack_input):
    '''This function receives a NUPACK input file and calculates the structure with the minimum free energy.
    It is important to remember that NUPACK input files must be named with the ".in" extension but this extension
    should not be included when calling NUPACK. For example: when executing NUPACK with a file called "example.in",
    we only need to say "example".
    '''
    command = 'mfe -material rna1995 -T 37 -multi ' + nupack_input
    print(command)
    a = os.system(command)
    return(a)

def parse_nupack(nupack_output):
    '''This function parses the output by NUPACK and returns a dictionary with the summary of the results,
    that is, the sequence, the free energy, the secondary structure, and the list of paired bases.
    '''
    handle = open(nupack_output, 'r')
    out_dict = OrderedDict()
    
    # This boolean will help us know if we are looking at the results
    # bool_results = False
    
    # This counter will help me load the data
    counter = 0
    
    paired_bases = []
    
    for line in handle:
        # Retrieve the sequence
        if 'Sequence:' in line:
            sequence = line[13:].strip()
        # Look for the separator that will start giving us the numbers we want
        elif line.startswith('% %'):
            # bool_results = True
            counter = 1
        elif counter == 1:
            # Then I know this is the length of the sequence.
            counter = counter + 1
        elif counter == 2:
            # Then I know this is the free energy
            energy = float(line)
            counter = counter + 1
        elif counter == 3:
            # This is the secondary structure
            structure = line.strip()
            counter = counter + 1
        elif counter > 3:
            # Then I am looking at the table of paired bases
            paired_bases.append(tuple(line.strip().split('\t')))
            
    handle.close()
    # Return the dictionary
    out_dict['sequence'] = sequence
    out_dict['energy'] = energy
    out_dict['structure'] = structure
    out_dict['paired_bases'] = paired_bases
    
    return(out_dict)

# The zip function goes through elements in matching positions of iterables until one of them runs out.
def hamming_distance(s1, s2):
    # assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
def make_blast_database (reference_genome_file, tag):
  '''
  Input: reference_genome_file the path to the input fasta file, the path to the output files which includes the base name for the files 
  output: creates blast databases in tag folder
  '''
  database_command = 'makeblastdb -in ' + reference_genome_file + ' -input_type fasta -dbtype nucl -out ' + reference_genome_file[0:-6]
  print(database_command)
  os.system(database_command)
 
def count_genome_matches(toeholds_fasta, working_directory, toehold_list, reference_genome_file, tag, new_df, pct_ident, evalue):
    '''This function will align toeholds to a reference genome and add the number of counts to a data frame'''
    #prepare Blast+ database 
    make_blast_database (reference_genome_file, tag)

    # Use the fasta file to align to the reference genome
    part1 = 'blastn -db ' + reference_genome_file[0:-6]
    part2 = ' -outfmt 6 -task blastn-short'
    part3 = ' -query ' + toeholds_fasta
    part4 = ' -perc_identity ' + pct_ident
    part5 = ' -evalue ' + evalue
    part6 = ' > ' + os.path.join(working_directory, tag + '_toeholds_alignment.aln')
    #print(part1 + part2 + part3 + part4 + part5 + part6)
    os.system(part1 + part2 + part3 + part4 + part5 + part6)
    #print('os.system(part1 + part2 + part3 + part4 + part5 + part6)')

    # Add the number of matches to the table
    handle = open(os.path.join(working_directory, tag + '_toeholds_alignment.aln'), 'r')
    reader = csv.reader(handle, delimiter = '\t')
    
    match_dictionary = OrderedDict()
    
    for line in reader:
        # Get the index of this toehold
        index = int(line[0][8:])
        
        # I could add a filter based on the alignment quality here
        
        if match_dictionary.get(index, -1) == -1:
            # Add this index to the dictionary
            match_dictionary[index] = 1
        else:
            # Add one match to this index
            match_dictionary[index] = match_dictionary[index] + 1
    
    handle.close()
    
    # Loop through the entries in the dataframe and add the number of matches
    column_name = 'Genome_matches_' + tag
    match_column = []
    
    # Loop through the dataframe to add the number of matches for each toehold
    for df_index, row in new_df.iterrows():
        #index = row[7]
        index = row[15]
        match_column.append(match_dictionary.get(int(index), 0))
    
    # Add the new column to the data frame
    new_df[column_name] = match_column
    
    return(new_df)
"""
def generate_toehold(trigger, mol_type, reporter, output_folder, length_unpaired):
    # Create the output directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Get the recognition sequence for our toehold
    recognition_sequence = get_switch_recognition_seq(trigger, mol_type, length_unpaired)
    
    # Design the starting sequence for the switch
    toehold = get_full_switch(recognition_sequence, reporter, length_unpaired)
    
    if toehold != 'Stop':
    
        # Use nupack to get the secondary structure for our switch
        first_toehold_nupack_input = os.path.join(output_folder, 'switch1_python.in')
        prepare_nupack_input(toehold, first_toehold_nupack_input)
        nupack_mfe(first_toehold_nupack_input[0:-3])
        
        # Read the secondary structure for the toehold we just generated
        first_toehold_mfe = os.path.join(output_folder, 'switch1_python.mfe')
        first_toehold_dict = parse_nupack(first_toehold_mfe)
"""
def generate_toehold(trigger, mol_type, reporter, output_folder, length_unpaired):
    # Create the output directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Get the recognition sequence for our toehold
    recognition_sequence = get_switch_recognition_seq(trigger, mol_type, length_unpaired)
    
    # Design the starting sequence for the switch
    toehold = get_full_switch(recognition_sequence, reporter, length_unpaired)
    
    if toehold != 'Stop':
    
        # Use nupack to get the secondary structure for our switch
        first_toehold_nupack_input = os.path.join(output_folder, 'switch1_python.in')
        prepare_nupack_input(toehold, first_toehold_nupack_input)
        nupack_mfe(first_toehold_nupack_input[0:-3])
        
        # Read the secondary structure for the toehold we just generated
        first_toehold_mfe = os.path.join(output_folder, 'switch1_python.mfe')
        first_toehold_dict = parse_nupack(first_toehold_mfe)

        # Use nupack to get the secondary structure for our RBSlinker region
        first_RBS_LINKER_nupack_input = os.path.join(output_folder, 'RBS_LINKER.in')
        RBS_LINKER = toehold[3+length_unpaired + length_paired-1: len(toehold)-1-len(reporter)] #pass length paired --> add - length reporter
        print (RBS_LINKER)
        prepare_nupack_input(RBS_LINKER, first_RBS_LINKER_nupack_input)
        nupack_mfe(first_RBS_LINKER_nupack_input[0:-3])
        
        # Read the secondary structure for the toehold we just generated
        first_RBS_LINKER_mfe = os.path.join(output_folder, 'RBS_LINKER.mfe')
        first_RBS_LINKER_dict = parse_nupack(first_RBS_LINKER_mfe)