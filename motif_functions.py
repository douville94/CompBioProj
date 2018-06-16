#!/usr/bin/python

######################################################
# Expectation-Maximization Algorithm for Finding Motifs
# Customizable # of sequences and motif lengths
# Spring 2018
# Trenton Langer
#######################################################

import random
import math
import time
import sys
import re

##################################################################
# generateSeqsFromFasta -   takes a list of FASTA organized sequences 
#                           and store the sequences in an array
# filename: name of the file to read/convert
#
# return: List of strings with amino acid characters
##################################################################
def generateSeqsFromFasta(filename):
    # read in file
    file = open(filename)

    # read in entire file
    allLines = file.readlines()
    seqs = []
    count = -1
  
    for i in range(len(allLines)):
        if(allLines[i][0] == '>'):
            seqs.append("")
            count += 1
        else:
            seqs[count] += allLines[i]
    
    for i in range(len(seqs)):
        seqs[i] = seqs[i].replace("\r", "")
        seqs[i] = seqs[i].replace("\n", "")
        
    # close file
    file.close()
    
    return seqs
    
##################################################################
# generateSeqsWithCutoffs -     takes a file of FASTA organized 
#                               sequences with cutoff points indicated 
#                               by a (##) at the end of the ">..." line 
#                               and returns the sequences from the FASTA 
#                               file cut to that length
# filename: name of file to read/convert
# returnType: string to indicate whether seq or its name will be returned
#
# return: list of sequences or list of names depending on parameter
##################################################################
def generateSeqsWithCutoffs(filename, returnType):
    # read in file
    file = open(filename)

    # read in entire file
    allLines = file.readlines()
    
    seqs = []
    seqNames = []
    count = -1
    seqCuts = []
  
    for i in range(len(allLines)):
        if(allLines[i][0] == '>'):
            seqs.append("")
            seqNames.append(allLines[i])
            count += 1
            m = re.search('>.*\(([0-9]*)\)', allLines[i])
            #print(m.group(1))
            seqCuts.append(m.group(1))
        else:
            seqs[count] += allLines[i]
    
    for i in range(len(seqs)):
        seqs[i] = seqs[i][0:int(seqCuts[i])]
        seqs[i] = seqs[i].replace("\r", "")
        seqs[i] = seqs[i].replace("\n", "")
        
    # close file
    file.close()
    
    if(returnType == "Sequences"):
        return seqs
    elif(returnType == "Names"):
        return seqNames
    else:
        print("Error: Invalid return type\n")
        return "Error"

##########################################################
# getBestMatchInAmino -     Suppose *sequence* is an Amino Acid sequence, and suppose
#                           profile is organised such as: profile[position in motif][amino acid], where
#                           profile[0][0] would be 'F1', profile[0][1] would be 'L1',
#                           profile[1][1] would be 'L2', etc...
#                           Order of Amino Acids is "FLIMVSPTAYHQNKDECWRG"
# sequence: amino acid string to check
# profile: motif profile to check against
# motifLength: length of desired motif
#
# return: string of length MotifLength that best matches profile
###########################################################
def getBestMatchInAmino(sequence, profile, motifLength):
    seq = sequence
    bestMatch = 0
    bestMatchIndex = 0
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']

    for i in range(len(sequence)-motifLength+1):
        tempMatch = 1
        
        for j in range(motifLength):
            for k in range(len(profile[j])):
                if seq[i+j] == aminoAcids[k]:
                    tempMatch = tempMatch*profile[j][k]
            
        if tempMatch > bestMatch:
            bestMatch = tempMatch
            bestMatchIndex = i
    
    return seq[bestMatchIndex:bestMatchIndex+motifLength]
    
###########################################################
# getFrequencyOfAminoAtPosition -   Calculates and returns the frequency of the given
#                                   amino acid, aa, at given position in list of sequences
#                                   If frequency is 0, return .001 (so info content calculation is
#                                   well-defined)
# aa: amino acid to search for (Char)
# position: indecx of sequence to check at
# sequences: list of amino acid strings to make profile from
#
# return: float percentage of aa frequency at index given
###########################################################
def getFrequencyOfAminoAtPosition(aa, position, sequences):
    count = 0.0
    numSeqs = len(sequences)
    
    for i in range(numSeqs):
        if sequences[i][position] == aa:
            count += 1
        
    rtnVal = count/numSeqs
    if rtnVal == 0:
        rtnVal += .001
    return rtnVal
    
###########################################################
# calcInfoContentAmino -    Calculates and returns the information content of
#                           a motif profile, given the values in the profile
#                           Assumes background frequency of 100/20 = 5% for each nucleotide
# profile: motif profile 
#
# return: float value of calculated information content
###########################################################
def calcInfoContentAmino(profile):
    # To calculate log base 2, use math.log(x, 2)
    sum = 0.0
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    for i in range(len(profile)):    #uses length of F row of profile table
        for j in range(len(aminoAcids)):
            profileVal = profile[i][j]
            sum += profileVal*math.log(profileVal/0.05, 2)
    
    return sum
    
#########################################################
# printMotifAmino - 
#                   Prints the specified motif model to the screen.
#                                     1       2       3       4       ...
#                             F       F1      F2      F3      F4      ...
#                             L       L1      L2      L3      L4      ...
#                             I       I1      I2      I3      I4      ...
#                             ...     ...     ...     ...     ...     ...
# profile: motifProfile
# motifLength: length of motif
#
# return: N/A
#########################################################
def printMotifAmino(profile, motifLength):
    line0 = ""
    lines = []
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    for i in range(len(aminoAcids)):
        lines.append(aminoAcids[i])

    for i in range(motifLength):
        line0 = line0 + "\t" + str(i)
        for j in range(len(lines)):
            lines[j] += "\t" + str("{0:.3f}".format(profile[i][j]))
    
    print(line0)
    for i in range(len(lines)):
        print(lines[i])
    print()
            
##########################################################
# findBestMotifAmino -  Uses expectation-maximization algorithm for finding
#                       the best motif in the any number of aa sequences
# sequences: list of sequences to find a motif in
# motifLength: desired motif length
# 
# return: motif profile generated from sequences
##########################################################
def findMotifAmino(sequences, motifLength):

    numSeqs = len(sequences)
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    
    # remember the motif instances from the previous iteration so we know when algorithm converges
    old_instance = []
    randomStart = []
    instance = []
    for i in range(numSeqs):
        old_instance.append("")
        randomStart.append(random.randint(0, len(sequences[i])-motifLength))
        #print("randStart = " + str(randomStart[i]))
        instance.append(sequences[i][randomStart[i]:randomStart[i]+motifLength])
    
    # repeat two steps of EM until convergence.
    convergence = False
    while(not convergence):    
        #temp seq list to test getfrequency
        seqs = []
        for i in range(numSeqs):
            seqs.append(instance[i])
        
        profile = []
        for i in range(motifLength):
            profile.append([])
            for j in range(len(aminoAcids)):
                profile[i].append(getFrequencyOfAminoAtPosition(aminoAcids[j], i, seqs))

        # print the motif model to the screen (comment out when running findBestMotif)
        #printMotifAmino(profile, motifLength)
        
        convergence = True
        for i in range(numSeqs):
            old_instance[i] = instance[i]   # re-assign old instances as the current instances to determine convergence
            instance[i] = getBestMatchInAmino(sequences[i], profile, motifLength)    # find best match in each sequence (step 2 of EM) 
            if(old_instance[i] != instance[i]):
                convergence = False
        
    return profile
    
##################################################################
# generateConsensusAmino -  Scans through profile and returns list of consensus sequences
# profile: motif profile
#
# return: list of potential consensus sequences
##################################################################
def generateConsensusAmino(profile):
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    consensusSeqs = []
    consensusSeqs.append("")
    for i in range(len(profile)):
        maxVal = max(profile[i])
        count = 0
        oldNumSeqs = 1
        for j in range(len(profile[i])):
            if profile[i][j] == maxVal:
                if count == 0:
                    for k in range(len(consensusSeqs)):
                        consensusSeqs[k] += aminoAcids[j]
                    count += 1
                    oldNumSeqs = len(consensusSeqs)
                else:
                    for k in range(oldNumSeqs):
                        cutoff = len(consensusSeqs[k]) - 1
                        consensusSeqs.append(consensusSeqs[k][0:cutoff] + aminoAcids[j])
                        count += 1
                    
        
    #print(consensusSeqs)
    #print(len(consensusSeqs))
    return consensusSeqs      
    
##################################################################
# consensusMatch -  Calculates the best match in list of consensus seqs and returns 
#                   information about best matches
# sequences: sequences to check for a best match
# consensusSeqs: list of generated possible consensus sequences
# returnInfo: decides what gets returned
# >>> "Score": scores
# >>> "Sequence": sequences
#
# return: list of scores/sequences representing best match info
##################################################################
def consensusMatch(sequences, consensusSeqs, returnInfo):
    matchScores = []
    matchSeqs = []
    for i in range(len(sequences)):
        bestScore = -10000
        bestConSeq = 0
        for j in range(len(consensusSeqs)):
            score = blosumScore(sequences[i], consensusSeqs[j])
            if score > bestScore:
                bestScore = score
                bestConSeq = j
                #print("Matching: " + sequences[i] + " with " + consensusSeqs[j])
        matchScores.append(bestScore)
        matchSeqs.append(consensusSeqs[bestConSeq])
    #print(matchScores)
    
    if(returnInfo == "Score"):
        return matchScores
    elif(returnInfo == "Sequence"):
        return matchSeqs
    else:
        print("ERROR: consensusMatch - invalid returnInfo parameter\n")
        return "ERROR"

##################################################################
# blosumScore -     Blosum Compare two AA sequences
# seq1: amino acid sequence 1
# seq2: amino acid sequence 2
#
# return: blosum62 score between sequences
##################################################################
def blosumScore(seq1, seq2):
    sum = 0
    mat = readMatrix("BLOSUM62.txt")
    
    #If seqs are not equal length, return 0
    if(len(seq1) != len(seq2)):
        return 0
    
    # note: mat is a list of 2 items: the 2D score array and
    # the 1D legend of amino acid characters
    # to get the index of "A" in the legend, you can write
    # mat[1].index("A")
    for x in range(len(seq1)):
        let1 = seq1[x]
        let2 = seq2[x]
        score = mat[0][mat[1].index(let1)][mat[1].index(let2)]
        sum = sum + score
    
    return sum

#############################################################
# readMatrix -  read data from file and store into score matrix
#               returns score matrix as 2D array and the amino acid table order
#               (completed for us)
# inputFile: the blosum/pam matrix to read in
# 
# return: 2D array representing matrix from txt file
#############################################################
def readMatrix(inputFile):
    # create 2D scoring matrix 24 x 24 (20 AAs, plus non-determinates)
    score = [[0 for x in range(24)] for x in range(24)]
    
    f = open(inputFile)

    # remove top six lines of comments plus row of AAs from BLOSUM 62 file
    for i in range(6):
        f.readline()

    # create array for the order of amino acids
    legend = f.readline().split()
    
    # assign values in score
    row = 0
    for line in f:
        # parse the line
        items = line.split()  #splits on whitespace
        # first value in item is the amino acid letter, so we want to
        # ignore it by starting index at 1
        for i in range(1,len(items)):
            score[row][i-1] = int(items[i])
        row+=1
    # return both the 2D table of scores and the 1D legend of the AA order
    return [score, legend]

##################################################################
# findBestMotifAmino -  Runs findMotifAmino "repetition" number of times and returns the motif
#                       with the highest information content
#                       When calling this function, be sure to comment out
#                       printing in the findMotif function
# sequences: list of aa sequences
# names: list of names corresponding to sequences
# motifLength: desired length of motif
# repetition: number of times to repeat running findMotifAmino
# seqPrinting: determine whether or not to print full sequence to console
#
# return: list of best matching subsequence for each sequence, their match score, and the consensus sequence it matches
##################################################################
def findBestMotifAmino(sequences, names, motifLength, repetition, seqPrinting):
    bestInfoContent = -10000000
    
    bestProfile = []
    instance = []
    
    totalTime = 0
    last10Time = 0
    for i in range(repetition):
        time1 = time.clock()
        m = findMotifAmino(sequences, motifLength)
        infoContent = calcInfoContentAmino(m)
        #print(infoContent)

        # keep best motif found so far
        if infoContent > bestInfoContent:
                bestProfile = m
                bestInfoContent = infoContent
        time2 = time.clock()
        totalTime += (time2-time1)
        last10Time += (time2-time1)
        if(i % 10 == 0 and True):
            #print(str(i) + "\t" + str(last10Time))
            progress = float(i*100 / repetition)
            etr = last10Time*(repetition/10)*(repetition-i)/(repetition)
            if i != 0:
                sys.stdout.write("Run Progress: %f%%, Estimated Time Remaining: %f seconds \r" % (progress, etr))
                sys.stdout.flush()
            else:
                sys.stdout.write("Run Progress: %f%% \r" % progress)
                sys.stdout.flush()
            last10Time = 0
        
	
    # display information about best motif
    m = bestProfile
    
    # print best motif
    printMotifAmino(m, motifLength)
    
    #generateConsensusAmino sequence
    consensusSeqs = generateConsensusAmino(m)
    
    # find best match in each sequence   
    for i in range(len(sequences)):
        instance.append(getBestMatchInAmino(sequences[i], m, motifLength))
        if(seqPrinting):
            print(sequences[i] + "\t" + instance[i])

    
    #get best consensus match score for each instance
    matchScores = consensusMatch(instance, consensusSeqs, "Score")
    matchSeqs = consensusMatch(instance, consensusSeqs, "Sequence")
    
    # print runtime
    print("Total Execution Time: " + "{0:.2f}".format(totalTime) + " seconds")
    
    # print info content
    print ("\nInformation content: " + str(bestInfoContent)+ "\n")
        
    for i in range(len(instance)):
        print(names[i] + "\t" + instance[i] + ", match score: " + str(matchScores[i]) + "\t(Consensus Sequence: " + str(matchSeqs[i]) + ")")
                
    saveToFile(m, names, instance, matchSeqs, matchScores, motifLength)
    return instance
  
##################################################################
# saveToFile -   gets run in findBestMotifAmino, saves console ouptput to 
#                text file for bio faculty
# profile: generated motif profile
# names: list of sequence names
# instance: list of best match instances
# conSensusSeq: list of matched consensus sequence
# matchScores: list of score of matched seq
# motifLength: length of motif
#
# No Return
# 
# Output File: "NT_Motifs_Results.txt" - saves version of console output for easy reading by bio faculty
##################################################################  
def saveToFile(profile, names, instance, conSensusSeq, matchScores, motifLength):
    file = open("NT_Motifs_Results.txt", "w")
    
    file.write("Output Results of Motif Finding Algorithm\n")
    file.write("->Motif Length = " + str(motifLength) + "\n")
    file.write("->NT sections with a sequence length of less than the Motif Length are not displayed\n\n")
    
    file.write("Generated Motif Profile:\n")
    line0 = ""
    lines = []
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    for i in range(len(aminoAcids)):
        lines.append(aminoAcids[i])

    for i in range(motifLength):
        line0 = line0 + "\t\t" + str(i)
        for j in range(len(lines)):
            lines[j] += "\t" + str("{0:.3f}".format(profile[i][j]))
    
    file.write(line0 + "\n")
    for i in range(len(lines)):
        file.write(lines[i] + "\n")
        
    file.write("\nInformation Content: " + str(calcInfoContentAmino(profile)) + "\n")
    file.write("\nBest Matching Motif in each sequence:\n")
    for i in range(len(names)):
        formattedName = re.search(r">(.*?)\(.*\)", names[i]).group(1)
        file.write("\t" + formattedName + "\n")
        file.write("\t-->Best Matching Subsequence: " + instance[i] + "\n")
        file.write("\t-->BLOSUM62 Match Score: " + str(matchScores[i]) + "\n")
        file.write("\t-->Consensus Sequence Matched: " + conSensusSeq[i] + "\n\n")
        
    file.close()
    

##################################################################
# findBestMotif -   Runs findMotif "repetition" number of times and returns the motif
#                   with the highest information content
#                   When calling this function, be sure to comment out
#                   printing in the findMotif function
# sequences: list of sequences
# motifLength: desired motiflength
# repetition: number of runs to perform EM on to get best result
# amino: boolean for whether giving AA or nucleotide sequences
# seqPrinting: whether or not full sequence is printed
#
# return: list of best matching indices
##################################################################
def findBestMotif(sequences, motifLength, repetition, amino, seqPrinting):
    bestInfoContent = -10000000
    
    bestProfile = []
    instance = []
    
    totalTime = 0
    last10Time = 0
    for i in range(repetition):
        time1 = time.clock()
        m = findMotif(sequences, motifLength, amino)
        infoContent = calcInfoContent(m)
        #print(infoContent)

        # keep best motif found so far
        if infoContent > bestInfoContent:
                bestProfile = m
                bestInfoContent = infoContent
        time2 = time.clock()
        totalTime += (time2-time1)
        last10Time += (time2-time1)
        if(i % 10 == 0 and True):
            print(str(i) + "\t" + str(last10Time))
            last10Time = 0
	
    # display information about best motif
    m = bestProfile
    
    # print best motif
    if(seqPrinting):
        printMotif(m, motifLength)

    # find best match in each sequence   
    for i in range(len(sequences)):
        instance.append(getBestMatchInSequence(sequences[i], m, motifLength))
        if(seqPrinting):
            print(sequences[i] + "\t" + instance[i])

    # print runtime
    print("Total Execution Time: " + "{0:.2f}".format(totalTime) + " seconds")
    
    # print info content
    print ("\nInformation content: " + str(bestInfoContent)+ "\n")
        
    for i in range(len(instance)):
        print(str(i) + "\t" + instance[i] + "\t" + base2amino(instance[i]))
        
    return instance

##################################################################
# amino2bases -     converts amino acids to nucleotide bases using 
#                   only one of the possible nucleotide sequences, not accounting
#                   for amino acids made by multiple base sequences
# aaSeq: an amino acid sequence
#
# return: a nucleotide sequence
##################################################################
def amino2bases(aaSeq):
    nucSeq = ""
    for i in range(len(aaSeq)):
        if aaSeq[i] == 'F':
            nucSeq += 'UUU'
        elif aaSeq[i] == 'L':
            nucSeq += 'UUA'
        elif aaSeq[i] == 'I':
            nucSeq += 'AUU'
        elif aaSeq[i] == 'M':
            nucSeq += 'AUG'
        elif aaSeq[i] == 'V':
            nucSeq += 'GUU'
        elif aaSeq[i] == 'S':
            nucSeq += 'UCU'
        elif aaSeq[i] == 'P':
            nucSeq += 'CCU'
        elif aaSeq[i] == 'T':
            nucSeq += 'ACU'
        elif aaSeq[i] == 'A':
            nucSeq += 'GCU'
        elif aaSeq[i] == 'Y':
            nucSeq += 'UAU'
        elif aaSeq[i] == 'H':
            nucSeq += 'CAU'
        elif aaSeq[i] == 'Q':
            nucSeq += 'CAA'
        elif aaSeq[i] == 'N':
            nucSeq += 'AAU'
        elif aaSeq[i] == 'K':
            nucSeq += 'AAA'
        elif aaSeq[i] == 'D':
            nucSeq += 'GAU'
        elif aaSeq[i] == 'E':
            nucSeq += 'GAA'
        elif aaSeq[i] == 'C':
            nucSeq += 'UGU'
        elif aaSeq[i] == 'W':
            nucSeq += 'UGG'
        elif aaSeq[i] == 'R':
            nucSeq += 'AGA'
        elif aaSeq[i] == 'G':
            nucSeq += 'GGU'
    return nucSeq
    
##################################################################
# bases2amino - complement of amino2bases to reverse the processing
# nucSeq: nucleotide sequence
#
# return: amino acid sequence
##################################################################
def base2amino(nucSeq):
    aaSeq = ""
    length = len(nucSeq)
    i = 0
    while (i + 3 <= length):
        if nucSeq[i:i+3] == 'UUU':
            aaSeq += 'F'
        elif nucSeq[i:i+3] == 'UUA':
            aaSeq += 'L'
        elif nucSeq[i:i+3] == 'AUU':
            aaSeq += 'I'
        elif nucSeq[i:i+3] == 'AUG':
            aaSeq += 'M'
        elif nucSeq[i:i+3] == 'GUU':
            aaSeq += 'V'
        elif nucSeq[i:i+3] == 'UCU':
            aaSeq += 'S'
        elif nucSeq[i:i+3] == 'CCU':
            aaSeq += 'P'
        elif nucSeq[i:i+3] == 'ACU':
            aaSeq += 'T'
        elif nucSeq[i:i+3] == 'GCU':
            aaSeq += 'A'
        elif nucSeq[i:i+3] == 'UAU':
            aaSeq += 'Y'
        elif nucSeq[i:i+3] == 'CAU':
            aaSeq += 'H'
        elif nucSeq[i:i+3] == 'CAA':
            aaSeq += 'Q'
        elif nucSeq[i:i+3] == 'AAU':
            aaSeq += 'N'
        elif nucSeq[i:i+3] == 'AAA':
            aaSeq += 'K'
        elif nucSeq[i:i+3] == 'GAU':
            aaSeq += 'D'
        elif nucSeq[i:i+3] == 'GAA':
            aaSeq += 'E'
        elif nucSeq[i:i+3] == 'UGU':
            aaSeq += 'C'
        elif nucSeq[i:i+3] == 'UGG':
            aaSeq += 'W'
        elif nucSeq[i:i+3] == 'AGA':
            aaSeq += 'R'
        elif nucSeq[i:i+3] == 'GGU':
            aaSeq += 'G'
        i += 3
    return aaSeq

##################################################################################################################
# Testing
##################################################################################################################

######################################
# test1 - makes sure generateSeqsFromFasta works by printing results to the console 
# No Parameters
# No Return
######################################
def test1():
    seqs = generateSeqsFromFasta("ceHaspins_CloseRelatives_ExpressionList.txt")
    for i in range(len(seqs)):
        print(seqs[i])
        
#test1()
    
######################################
# testConsensus - mock up a profile and check what consensus sequences are generated by generateConsensusAmino()
# No Parameters
# No Return
######################################
def testConsensus():
    profile = []
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    for i in range(3):
        profile.append([])
        for j in range(len(aminoAcids)):
            if(j < 3):
                profile[i].append(0.25)
            else:
                profile[i].append(0.1)
            
    profile[2][5] = 0.3
    profile[1][5] = 0.3
    
    printMotifAmino(profile, len(profile))
    consensusSeqs = generateConsensusAmino(profile)
    print("Con Seqs: " + str(consensusSeqs))
    
    testSeqs = ["FLLLLL", "FSS", "AYHNSDNNKSD", "TAYQLSS"]
    
    # find best match in each sequence   
    instance = []
    for i in range(len(testSeqs)):
        instance.append(getBestMatchInAmino(testSeqs[i], profile, 3))
        print(instance[i])
    
    #get best consensus match score for each instance
    matchScores = consensusMatch(instance, consensusSeqs, "Score")
    matchSeqs = consensusMatch(instance, consensusSeqs, "Sequence")
    
    for i in range(len(instance)):
        print(str(i) + "\t" + instance[i] + ", match score: " + str(matchScores[i]) + "(" + str(matchSeqs[i]) + ")")
         
#testConsensus()

######################################
# testSeqCuts - tests generateSeqsWithCutoffs() by printing result
# No Parameters
# No Return
######################################
def testSeqCuts():
    aaSeqs = generateSeqsWithCutoffs("new_amino_acids.txt")
    print(aaSeqs)
    
#testSeqCuts()

######################################
# testAnalyzeHaspins - run algorithm without using generated NT cutoff point
# No Parameters
# No Return
######################################
def testAnalyzeHaspins(filename, motifLength, numberRuns):
    aaSeqs = generateSeqsFromFasta(filename)
    printFullSeqs = False
    instances = findBestMotifAmino(aaSeqs, motifLength, numberRuns, printFullSeqs)

#testAnalyzeHaspins("ceHaspins_CloseRelatives_ExpressionList.txt", 5, 10)

##################################################################################################################
# Run Functions Here
##################################################################################################################
######################################
# analyzeWithCutoffs -  reads an input file with cutoff information, extracts the sequences,
#                       runs findBestMotifAmino() and saves results to a .txt
# filename: name of file with amino acids and cutoff information
# motifLength: length of desired motif
# numberRuns: number of repetitions to put into findBestMotifAmino
#
# No Return 
######################################
def analyzeWithCutoffs(filename, motifLength, numberRuns):
    aaSeqs = generateSeqsWithCutoffs(filename, "Sequences")
    aaSeqsNames = generateSeqsWithCutoffs(filename, "Names")
    relevantSeqs = []
    relevantNames = []
    for i in range(len(aaSeqs)):
        if len(aaSeqs[i]) > motifLength:
            relevantSeqs.append(aaSeqs[i])
            relevantNames.append(aaSeqsNames[i])
    printFullSeqs = False
    instances = findBestMotifAmino(relevantSeqs, relevantNames, motifLength, numberRuns, printFullSeqs)
    
analyzeWithCutoffs("new_amino_acids.txt", 10, 10000)   

