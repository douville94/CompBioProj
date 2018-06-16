
########################################################
# convertFileToSequence - takes a FASTA file and returns the
# DNA sequence as a string
#######################################################
def convertFileToSequence(filename):
    # read in file
    file = open(filename)
    
    # read in first line
    header = file.readline()
    if (header[0] == '>'):
        print("in FASTA format")
    else:
        print("invalid format")
        file.close()
        # error code of -1 (so functions that call this can check)
        return -1
    
    # read in rest of file
    sequence = file.read()
    
    # close file
    file.close()

    # remove all return and newline characters
    sequence = sequence.replace("\r", "")
    sequence = sequence.replace("\n", "")
    return sequence

########################################################
# translate -- given the input file, translate each DNA
# codon to the corresponding amino acid sequence and writing
# it out to the given output file
#######################################################
def translate(inputfile, outputfile):
    # convert the input file to a string
    string = convertFileToSequence(inputfile)
    
    # open the file to write to
    out = open(outputfile, "w")
    out.write("> " + inputfile + " translated to amino acids")
    
    # counter to keep track of how many characters were written to the file
    count = 0
    
    # identifies the codon and finds the corresponding single letter amino acid
    for i in range(0, len(string), 3):
        if count % 60 == 0:
            out.write("\n")
        codon = string[i:i+3]
        if(codon == "TTT" or codon =="TTC" ):
            out.write("F")
            count += 1
        elif(codon == "TTA" or codon =="TTG" or codon =="CTT" or codon =="CTC" or codon == "CTA" or codon =="CTG"):
            out.write( "L")
            count += 1
        elif(codon == "ATT" or codon =="ATC" or codon =="ATA"):
            out.write( "I")
            count += 1
        elif(codon == "GTT" or codon =="GTC" or codon =="GTA"or codon =="GTG"):
            out.write( "V")
            count += 1
        elif(codon == "TTA" or codon =="TTG" or codon =="CTT"or codon =="CTC"or codon =="CTA" or codon =="CTG"):
            out.write( "L")
            count += 1
        elif(codon == "ATG" ):
            out.write("M")
            count += 1
        elif(codon == "TCT" or codon =="TCC" or codon =="TCA" or codon =="TCG"):
            out.write("S")
            count += 1
        elif(codon == "CCT" or codon =="CCC" or codon =="CCA" or codon =="CCG"):
            out.write("P")
            count += 1
        elif(codon == "ACT" or codon =="ACC" or codon =="ACA" or codon =="ACG"):
            out.write( "T" )
            count += 1
        elif(codon == "GCT" or codon =="GCC" or codon =="GCA" or codon =="GCG"):
            out.write( "A")
            count += 1
        elif(codon == "TAT" or codon =="TAC"):
            out.write( "Y")
            count += 1
        elif(codon == "TAA" or codon =="TAG" or codon =="TGA"):
            out.write( "X")
            count += 1
        elif(codon == "CAT" or codon =="CAC"):
            out.write( "H")
            count += 1
        elif(codon == "CAA" or codon =="CAG"):
            out.write( "Q")
            count += 1
        elif(codon == "AAT" or codon =="AAC" ):
            out.write( "N")
            count += 1
        elif(codon == "AAA" or codon =="AAG" ):
            out.write( "K")
            count += 1
        elif(codon == "GAT" or codon =="GAC"):
            out.write("D")
            count += 1
        elif(codon == "GAA" or codon =="GAG"):
            out.write("E")
            count += 1
        elif(codon == "TGT" or codon =="TGC"):
            out.write("C")
            count += 1
        elif(codon == "TGG"):
            out.write( "W" )
            count += 1
        elif(codon == "CGT" or codon =="CGC" or codon =="CGA" or codon =="CGG"):
            out.write( "R")
            count += 1
        elif(codon == "AGT" or codon =="AGC"):
            out.write( "S")
            count += 1
        elif(codon == "AGA" or codon =="AGG" ):
            out.write( "R")
            count += 1
        elif(codon == "GGT" or codon =="GGC" or codon =="GGA" or codon =="GGG"):
            out.write("G")
            count += 1

    # close the file
    out.close()


##################################################################
# generateSeqsFromFasta - takes a list of FASTA organized sequences
# and store the sequences in an array
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
# findIndex - finds the index of the beginning of the C-terminus and
# writes it to the end of the first line for each organism in the file
# 'new_amino_acids.txt'
##################################################################
def findIndex():
    translate("human_haspin_FASTA.txt", "human_haspin_amino_acid.txt")
    human_file = open("human_haspin_amino_acid.txt")

    human_amino_acid = convertFileToSequence("human_haspin_amino_acid.txt")
    amino_acid = generateSeqsFromFasta("amino_acids.txt")

    print()
    human = 0
    human_new = human_amino_acid.find("VCF", human)
    if human_new == 0:
        human = 0
    else:
        human = len(human_amino_acid)
    while human != human_new:
        if human_new < human:
            human = human_new
            if human_amino_acid.find("VCF", human) != -1:
                human_new = human_amino_acid.find("VCF", human)
    print("human_amino_acid", human)

    C01H6 = 0
    new_C01H6 = amino_acid[0].find("MST", C01H6)
    if new_C01H6 == 0:
        C01H6 = 0
    else:
        C01H6 = len(amino_acid[0])
    while C01H6 != new_C01H6:
        if new_C01H6 < C01H6:
            C01H6 = new_C01H6
            if amino_acid[0].find("MST", C01H6) != -1:
                new_C01H6 = amino_acid[0].find("MST", C01H6)
    print("C01H6.9_hasp1_conceptual", C01H6)

    Y18H1A = 0
    new_Y18H1A = amino_acid[1].find("WMR", Y18H1A)
    if new_Y18H1A == 0:
        Y18H1A = 0
    else:
        Y18H1A = len(amino_acid[1])
    while Y18H1A != new_Y18H1A:
        if new_Y18H1A < Y18H1A:
            Y18H1A = new_Y18H1A
            if amino_acid[1].find("WMR", Y18H1A) != -1:
                new_Y18H1A = amino_acid[1].find("WMR", Y18H1A)
    print("Y18H1A.10_hasp2_conceptual", Y18H1A)

    H12I13 = 0
    new_H12I13 = amino_acid[7].find("KCQ", H12I13)
    if new_H12I13 == 0:
        H12I13 = 0
    else:
        H12I13 = len(amino_acid[7])
    while H12I13 != new_H12I13:
        if new_H12I13 < H12I13:
            H12I13 = new_H12I13
            if amino_acid[7].find("KCQ", H12I13) != -1:
                new_H12I13 - amino_acid[7].find("KCQ", H12I13)
    print("Y18H1A.H12I13.1conceptual", H12I13)

    ZK177 = 0
    new_ZK177 = amino_acid[2].find("LCI", ZK177)
    if new_ZK177 == 0:
        ZK177 = 0
    else:
        ZK177 = len(amino_acid[2])
    while ZK177 != new_ZK177:
        if new_ZK177 < ZK177:
            ZK177 = new_ZK177
            if amino_acid[2].find("LCI", ZK177) != -1:
                new_ZK177 = amino_acid[2].find("LCI", ZK177)
    print("ZK177.2conceptual", ZK177)

    F59E12 = 0
    new_F59E12 = amino_acid[4].find("KCD", F59E12)
    if new_F59E12 == 0:
        F59E12 = 0
    else:
        F59E12 = len(amino_acid[4])
    while F59E12 != new_F59E12:
        if new_F59E12 < F59E12:
            F59E12 = new_F59E12
            if amino_acid[4].find("KCD", F59E12) != -1:
                new_F59E12 = amino_acid[4].find("KCD", F59E12)
    print("F59E12.6and15_combined", F59E12)

    C55C3 = 0
    new_C55C3 = amino_acid[9].find("MAW", C55C3)
    if new_F59E12 == 0:
        C55C3 = 0
    else:
        C55C3 = len(amino_acid[9])
    while C55C3 != new_C55C3:
        if new_C55C3 < C55C3:
            C55C3 = new_C55C3
            if amino_acid[9].find("MAW", C55C3) != -1:
                new_C55C3 = amino_acid[9].find("MAW", C55C3)
    print("C55C3.8conceptual", C55C3)

    Y73B6A = 0
    new_Y73B6A = amino_acid[5].find("QMM", Y73B6A)
    if new_Y73B6A == 0:
        Y73B6A = 0
    else:
        Y73B6A = len(amino_acid[5])
    while Y73B6A != new_Y73B6A:
        if new_Y73B6A < Y73B6A:
            Y73B6A = new_Y73B6A
            if amino_acid[5].find("QMM", Y73B6A) != -1:
                new_Y73B6A = amino_acid[5].find("QMM", Y73B6A)
    print("Y73B6A.1aconceptual", Y73B6A)

    C04G2 = 0
    new_C04G2 = amino_acid[3].find("MSP", C04G2)
    if new_C04G2 == 0:
        C04G2 = 0
    else:
        C04G2 = len(amino_acid[3])
    while C04G2 != new_C04G2:
        if new_C04G2 < C04G2:
            C04G2 = new_C04G2
            if amino_acid[3].find("MSP", C04G2) != -1:
                new_C04G2 = amino_acid[3].find("MSP", C04G2)
    print("C04G2.10aconceptual", C04G2)

    Y32G9A = 0
    new_Y32G9A = amino_acid[11].find("MEH", Y32G9A)
    if new_Y32G9A == 0:
        Y32G9A = 0
    else:
        Y32G9A = len(amino_acid[11])
    while Y32G9A != new_Y32G9A:
        if new_Y32G9A < Y32G9A:
            Y32G9A = new_Y32G9A
            if amino_acid[11].find("MEH", Y32G9A) != -1:
                new_Y32G9A = amino_acid[11].find("MEH", Y32G9A)
    print("Y32G9A.3conceptual", Y32G9A)

    C26E6 = 0
    new_C26E6 = amino_acid[10].find("MKG", C26E6)
    if new_C26E6 == 0:
        C26E6 = 0
    else:
        C26E6 = len(amino_acid[10])
    while C26E6 != new_C26E6:
        if new_C26E6 < C26E6:
            C26E6 = new_C26E6
            if amino_acid[10].find("MKG", C26E6) != -1:
                new_C26E6 = amino_acid[10].find("MKG", C26E6)
    print("C26E6.1conceptual", C26E6)

    Y48B6A = 0
    new_Y48B6A = amino_acid[6].find("KGK", Y48B6A)
    if new_Y48B6A == 0:
        Y48B6A = 0
    else:
        Y48B6A = len(amino_acid[6])
    while Y48B6A != new_Y48B6A:
        if new_Y48B6A < Y48B6A:
            Y48B6A = new_Y48B6A
            if amino_acid[6].find("KGK", Y48B6A) != -1:
                new_Y48B6A = amino_acid[6].find("KGK", Y48B6A)
    print("Y48B6A.10conceptual", Y48B6A)

    C50H2 = 0
    new_C50H2 = amino_acid[8].find("MTP", C50H2)
    if new_C50H2 == 0:
        C50H2 = 0
    else:
        C50H2 = len(amino_acid[8])
    while C50H2 != new_C50H2:
        if new_C50H2 < C50H2:
            C50H2 = new_C50H2
            if amino_acid[8].find("MTP", C50H2) != -1:
                new_C50H2 = amino_acid[8].find("MTP", C50H2)
    print("C50H2.7conceptual", C50H2)


    file = open("new_amino_acids.txt", "w")
    file2 = open("amino_acids.txt", "r")
    for line in file2:
        if line == ">C01H6.9_hasp1_conceptual translation\n":
            file.write(">C01H6.9_hasp1_conceptual translation (" + str(C01H6) + ")\n")
            file.write(amino_acid[0] + "\n\n")
        if line == ">Y18H1A.10_hasp2_conceptual translation\n":
            file.write(">Y18H1A.10_hasp2_conceptual translation (" + str(Y18H1A) + ")\n")
            file.write(amino_acid[1] + "\n\n")
        if line == ">ZK177.2conceptual translation\n":
            file.write(">ZK177.2conceptual translation (" + str(ZK177) + ")\n")
            file.write(amino_acid[2] + "\n\n")
        if line == ">C04G2.10aconceptual translation\n":
            file.write(">C04G2.10aconceptual translation (" + str(C04G2) + ")\n")
            file.write(amino_acid[3] + "\n\n")
        if line == ">F59E12.6and15_combined\n":
            file.write(">F59E12.6and15_combined (" + str(F59E12) + ")\n")
            file.write(amino_acid[4] + "\n\n")
        if line == ">Y73B6A.1aconceptual translation\n":
            file.write(">Y73B6A.1aconceptual translation (" + str(Y73B6A) + ")\n")
            file.write(amino_acid[5] + "\n\n")
        if line == ">Y48B6A.10conceptual translation\n":
            file.write(">Y48B6A.10conceptual translation (" + str(Y48B6A) + ")\n")
            file.write(amino_acid[6] + "\n\n")
        if line == ">H12I13.1conceptual translation\n":
            file.write(">H12I13.1conceptual translation (" + str(H12I13) + ")\n")
            file.write(amino_acid[7] + "\n\n")
        if line == ">C50H2.7conceptual translation\n":
            file.write(">C50H2.7conceptual translation (" + str(C50H2) + ")\n")
            file.write(amino_acid[8] + "\n\n")
        if line == ">C55C3.8conceptual translation\n":
            file.write(">C55C3.8conceptual translation (" + str(C55C3) + ")\n")
            file.write(amino_acid[9] + "\n\n")
        if line == ">C26E6.1conceptual translation\n":
            file.write(">C26E6.1conceptual translation (" + str(C26E6) + ")\n")
            file.write(amino_acid[10] + "\n\n")
        if line == ">Y32G9A.3conceptual translation\n":
            file.write(">Y32G9A.3conceptual translation (" + str(Y32G9A) + ")\n")
            file.write(amino_acid[11] + "\n\n")

    file.close
    file2.close


##################################################################
# testing
##################################################################
findIndex()














