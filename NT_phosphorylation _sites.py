"""
NT_phosphorolation_sites    Finds motifs and motif occurances containing
                            phosphorolation sites.

Author: Megan Nalani Chun
"""

import random
import math

"""
hertzStormo     Finds motifs based on a given singleton to begin with

parameters     
    aminoAcid - The amino acid the singleton has. This is a phosphorolation
                site (Either S: serine, Y:tyrosine, T:threoine)
    singleton - The base or starting point for hertz stormo
    singletonSeq - The index for the sequence that the singleton was selected from
    singletonI - The index where the amino acid is in the sequence with the singleton
    seqs - All sequences in the input file that could be matched with a motif
    headers - All headers for all the sequences

return
    bestSet - The set of motif matches
    bestHeaders - The set of headers where the motif matches were found
"""
def hertzStormo(aminoAcid, singleton, singletonSeq, singletonI, seqs, headers):
    size = len(singleton)                               #motif size
    bestIC = 0.0                                        #best information content
    bestSet = [singleton]                               
    bestHeaders = [headers[singletonSeq]]               
    bestLocations = [singletonI]                        #index locations for matches
    doneSeqs = [singletonSeq]                           #seqences already looked at
    mini = math.floor(float(len(bestSet[0]))/3.0)+1     #last check for matches
    if mini < 2:
        mini = 2

    #run hertz stormo algorithm, each iteration adds another match to the set until all sequences are considered
    for x in range(len(seqs)-1):
        addIC, addSet, addHeaders, addLocations, addSeq = iteration(seqs, headers, bestSet, doneSeqs, aminoAcid, size)
        same = 0
        for c in range(len(addSet)):
            if addSet[c] == bestSet[0][c]:
                same += 1
        if same >= mini:
            bestHeaders.append(addHeaders)
            bestIC = addIC
            bestLocations.append(addLocations)
            bestSet.append(addSet)
            doneSeqs.append(addSeq)

    return bestSet, bestHeaders

"""
iteration   Finds the next best motif match to add to the current
            set of best motifs based on the information content 
            when trying all possible N-mers from sequences that
            haven't been already looked at as long as they contain
            the original amino acid phosphorolation site. Helper
            for hertzStormo

parameters
    seqs - All sequences in the input file that could be matched with a motif
    headers - All headers for all the sequences
    currentSet - The current set of motif matches that will be part of the output
    doneSeqs - Indexes of the sequences that have already been considered
    aminoAcid - The amino acid the singleton has. This is a phosphorolation
                site (Either S: serine, Y:tyrosine, T:threoine) 
    size - The size of the prospective motifs to look at

return 
    bestIC - The best information content for this match
    bestNMer - The best match found
    bestHeader - The header for the sequence the match is from
    bestLocation - The index of the amino acid in the sequence
    bestSeq - The index for the sequence the match is from
"""
def iteration(seqs, headers, currentSet, doneSeqs, aminoAcid, size):    
    bestNMer = ""
    bestHeader = ""
    bestLocation = ""
    bestIC = 0.0
    bestSeq = 0
    
    for s in range(len(seqs)):
        #only search new sequences
        if s not in doneSeqs:
            i = 0                
            while i != -1:
                #only test match prospects with the amino acid
                i = seqs[s].find(aminoAcid, i+1)
                if i != -1:
                    #test all possible n lengths with the amino acid
                    for x in range(size):
                        #bounds check
                        if i-x >= 0 and i+size-x < len(seqs[s]):
                            #test match prospect only with the chosen set 
                            temp = []
                            for c in currentSet:
                                temp.append(c)
                            temp.append(seqs[s][i-x:i+size-x])

                            #get information content from profile to compare
                            tempProfile = createProfile(size, temp)
                            val = getIC(tempProfile)

                            #only remember the best match
                            if val > bestIC:
                                bestIC = val
                                bestNMer = seqs[s][i-x:i+size-x]
                                bestHeader = headers[s]                            
                                bestLocation = i-x
                                bestSeq = s
    return bestIC, bestNMer, bestHeader, bestLocation, bestSeq

"""
getIC   Calculates the information content (IC)

parameter
    profile - The profile to calculate the IC for

return 
    i - The IC
"""
def getIC(profile):
    #background = [.085, .095, .07, .01, .06, .095, .06, .02, .005, .04, .05, .06, .03, .03, .04, .06, .05, .02, .06, .06]    
    background = [.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]    
    i = 0
    for col in range(len(profile[0])):        
        for row in range(len(profile)):
            check = profile[row][col]
            if check > 0.0:
                i += check * math.log(check/background[row])
    return i
                    
"""
createProfile   Build the motif profile

parameters
    numCols - Number of columns for the profile
    seqs - The sequences to build the profile from

return
    profile - The motif profile for the given sequences
"""
def createProfile(numCols, seqs):
    profile = []  
    for col in range(numCols):            
        avgs = getAvg(seqs, col)
        profile.append(avgs)
    return profile

"""
getAvg  Calculates the average percent of occurances for each amino acid for. 
        every column in the motif profile. Helper for createProfile.

parameters
    seqs - The sequences to calculate the averages for
    col - The column of the motif profile to look at

return
    avgs - The averages for each amino acid occurance for
           the given column in a motif profile
"""
def getAvg(seqs, col):
    avgs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for s in seqs:
        aminoAcid = s[col]
        if aminoAcid == 'G':
            avgs[0] += 1
        elif aminoAcid == 'A':
            avgs[1] += 1
        elif aminoAcid == 'V':
            avgs[2] += 1
        elif aminoAcid == 'C':
            avgs[3] += 1
        elif aminoAcid == 'P':
            avgs[4] += 1
        elif aminoAcid == 'L':
            avgs[5] += 1
        elif aminoAcid == 'I':
            avgs[6] += 1
        elif aminoAcid == 'M':
            avgs[7] += 1
        elif aminoAcid == 'W':
            avgs[8] += 1
        elif aminoAcid == 'F':
            avgs[9] += 1
        elif aminoAcid == 'S':
            avgs[10] += 1
        elif aminoAcid == 'T':
            avgs[11] += 1
        elif aminoAcid == 'Y':
            avgs[12] += 1
        elif aminoAcid == 'N':
            avgs[13] += 1
        elif aminoAcid == 'Q':
            avgs[14] += 1
        elif aminoAcid == 'K':
            avgs[15] += 1
        elif aminoAcid == 'R':
            avgs[16] += 1
        elif aminoAcid == 'H':
            avgs[17] += 1
        elif aminoAcid == 'D':
            avgs[18] += 1
        elif aminoAcid == 'E':
            avgs[19] += 1
        else:
            print(aminoAcid, " is not accounted for")

    for i in range(20):
        avgs[i] = avgs[i]/len(seqs)
    return avgs

"""
findPhosSiteMotif   Finds a singleton for the hertz stormo algorithm. The
                    singleton must have a phosphorolation site.

parameteres
    sequences - List of sequences to search
    headers - List of headers for the sequences
    longest - Longest singleton that should be considered
    shortest - Shortest singleton that should be considered

return
    s[aminoAcid] - The amino acid phosphorolation site
    motif - The singleton that was chosen
    motifSeq - The index of the sequence the motif was found in
    aminoAcid - The index location in the sequence the amino acid is
    motifSeqs - The sequences that should be used for hertz stormo
    motifHeaders  - The headers for the motifSeqs
"""
def findPhosSiteMotif(sequences, headers, longest, shortest):    
    seqsLen = len(sequences)
    for i in range(seqsLen + 3):
        #randomly select a sequence to look at
        r = random.randint(0, seqsLen-1)
        s = sequences[r]     
        if len(s) >= shortest:
            for aminoAcid in range(len(s)):                
                #get singleton of proper length if phosphorolation site
                if s[aminoAcid] == 'S' or s[aminoAcid] == 'Y' or s[aminoAcid] == 'T':
                    motif = getMotif(s, s[aminoAcid], longest, shortest)
                    l = len(motif)
                    if l >= shortest:
                        #align sequences and headers to use for hertz stormo
                        motifSeqs = []
                        motifHeaders = []
                        motifSeq = 0
                        unaccounted = 0
                        for x in range(seqsLen):
                            if len(sequences[x]) >= l:                                                                                                                                                            
                                motifSeqs.append(sequences[x])
                                motifHeaders.append(headers[x])
                                if s == sequences[x]:                                    
                                    motifSeq = x-unaccounted                                
                            else:
                                unaccounted += 1
                        return s[aminoAcid], motif, motifSeq, aminoAcid, motifSeqs, motifHeaders    
    return None, None, None

"""
getMotif    Randomly selects a singleton. 
            Helper for findPhosSiteMotif

parameters
    seq - Sequence to look at
    aminoAcid - The index location in the sequence the amino acid is
    longest - Longest singleton that should be considered
    shortest - Shortest singleton that should be considered

return 
    seq[i-front:i+back] - A singleton
"""
def getMotif(seq, aminoAcid, l , s):
    i = seq.index(aminoAcid)
    front = 500
    back = 500
    while ((i-front) < 2) or ((i+back) >= len(seq)) or len(seq[i-front:i+back]) < s:
        front = random.randint(0,l)
        back = random.randint(0,l)
    return seq[i-front:i+back]

"""
getNumPhosphorolationSites  Print the number of phosphorolation sites of each type

parameter
    sequences - All sequences in input file
"""
def getNumPhosphorolationSites(sequences):
    sAA = []
    yAA = []
    tAA = []
    for s in sequences:
        sI = 0
        yI = 0
        tI = 0
        #get number of S amino acids in input file
        while sI != -1 and sI+1 < len(s):
            if sI == 0:
                sI = s.find('S', sI)
            else:
                sI = s.find('S', sI+1)
            if sI != -1:
                sAA.append(sI)

        #get number of Y amino acids in input file
        while yI != -1 and yI+1 < len(s):
            if yI == 0:
                yI = s.find('Y', yI)
            else:
                yI = s.find('Y', yI+1)
            if yI != -1:
                yAA.append(yI)

        #get number of T amino acids in input file
        while tI != -1 and tI+1 < len(s):
            if tI == 0:
                tI = s.find('T', tI)
            else:
                tI = s.find('T', tI+1)
            if tI != -1:
                tAA.append(tI)
    print("SAA: ", sAA)
    print("YAA: ", yAA)
    print("TAA: ", tAA)

#check if a line is a header
def isHeader(line):
    if (line[0] == '>'):
        return True
    else:
        return False

"""
convertFileToSequences  Converts a file in FASTA format containing multiple
                        sequences of amino acids into a list of the 
                        sequences and a list of the headers. File must have
                        the index for the N terminal half of the sequence
                        at the end of the header. Ex: >header (20)

parameter
    filename - The file to retreive data from

return
    seqs - The sequences in the file
    headers - The headers for the seqs    
"""
def convertFileToSequences(filename):
    seqs = []
    headers = []

    with open(filename) as f: 
        line = f.readline()        

        while line:
            #check for N-terminal index
            tempSeq = ""
            i = line.index('(')
            if i == -1:
                print("Sorry, file is in wrong format. Exiting now.")
                return
            #get index for N terminal
            NT_end = int(line[i+1:len(line)-2])
            headers.append(line[:len(line)-1])
            line = f.readline()
            if not line:
                break
            #get next sequence
            while not isHeader(line):
                tempSeq += line.strip()
                line = f.readline()
                if not line:
                    break
        
            #only keep N-terminal half of the sequence
            seqs.append(tempSeq[:NT_end])

    return seqs, headers

#get final motif for all the matches after running hertz stormo
def getFinal(matches):
    aa = ["G", "A", "V", "C", "P", "L", "I", "M", "W", 
        "F", "S", "T", "Y", "N", "Q", "K", "R", "H", "D", "E"]
    motif = ""
    profile = createProfile(len(matches[0]), matches)

    for col in range(len(profile)):
        largest = 0.0
        largestA = 0
        for aminoAcid in range(len(profile[0])):
            if profile[col][aminoAcid] > largest:
                largest = profile[col][aminoAcid]
                largestA = aminoAcid
        motif += aa[largestA] 
    return motif

#write header to output file
def setHeader(filename, inputfile):
    with open(filename, "w") as f:
        f.write("Output for NT phosphorolation sites from " + inputfile + "\n")

#write output to output file
def writeData(filename, sets, headers, motif):
    with open(filename, 'a') as f:
        f.write("\n----------- motif " + motif + "-----------\n")
        for i in range(len(sets)):
            f.write(headers[i] + "\n" + sets[i] + "\n")

#prompt user to provide parameters for program
inputfile = input("What data file would you like to use? ")
outputfile = "output_NT_phos_" + inputfile
longestMotifSize = int(input("What is the longest motif length you want? "))
shortestMotifSize = int(input("What is the shortest motif length you want? "))
numberOfMatches = int(input("What is the least number of occurances you are intested in for each motif? "))
while longestMotifSize < shortestMotifSize or longestMotifSize < 1 or shortestMotifSize < 1 or numberOfMatches < 1:
    print("\nSorry, the last three inputs don't make sense. Try again.\n")
    longestMotifSize = int(input("What is the longest motif length you want? "))
    shortestMotifSize = int(input("What is the shortest motif length you want? "))
    numberOfMatches = int(input("What is the least number of occurances you are intested in for each motif? "))

#prepare the output file
setHeader(outputfile, inputfile)

#extract the data
seqs, headers = convertFileToSequences(inputfile)

#Run hertx stormo with various singletons to get motifs and their matches
motifsSoFar = []
for x in range(100):
    #get a singleton
    aminoAcid, motif, motifSeq, motifI, seqs, headers = findPhosSiteMotif(seqs, headers, longestMotifSize, shortestMotifSize)

    #run hertz stormo and keep as output if in bounds
    if motif not in motifsSoFar:
        bestSet, bestHeaders = hertzStormo(aminoAcid, motif, motifSeq, motifI, seqs, headers)
        if len(bestSet) >= numberOfMatches:
            finalMotif = getFinal(bestSet)
            writeData(outputfile, bestSet, bestHeaders, finalMotif)
        motifsSoFar.append(motif)

print("\nDone finding motifs. To see your results, please open " + outputfile + "\n")