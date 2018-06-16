#IC and motifs
#!/usr/bin/python

#Luke Douville's part
#Look for motifs in the N terminus using algorithms that generate motifs.

#Locate the N terminus - it's at the beginning of the amino acid sequence
#Run algorithm to find motifs

import motif_functions

aaFile = open("amino_acids.txt", "r")
aaText = aaFile.read()
d1 = {}
def isolateAminoAcidSequences():#x, y):#aaString):
    global aaString
    aaString = ""
    for i in range(0, len(aaText)):
        if i + 1 < aaText.find(">", i, len(aaText)):
            aaString += aaText[i]#capture one amino acid sequence
        elif aaText[i] == "\n" and aaText[i + 1] == ">" and i < len(aaText):
            global d1
            d1[i] = aaString
            aaString = ""
        else:#quit capturing amino acid sequence when the end of the sequence is reached
            aaString = ""

    for i in range(0, len(aaText)):
        if i in d1:
            print(i, " ", end = ",")
    print("\n\n", d1[5324])
    print("\n\n\n")

isolateAminoAcidSequences()#Each amino acid sequence is now in d1
#print("d1 = ", d1)

#Convert amino acid sequences to RNA
d2 = {}
n = 1
for i in range(len(aaText)):
    if i in d1:
        d2[i] = motif_functions.amino2bases(d1[i])
        print("\nRNA sequence from amino acid sequence {} (in d1): ".format(n), motif_functions.amino2bases(d1[i]))
        n += 1

o = 1
for i in range(len(aaText)):
    if i in d2:
        print("\nRNA sequence from amino acid sequence {} (in d2): ".format(o), motif_functions.amino2bases(d2[i]))
        o += 1

#Look at the beginning of the DNA sequence for the best motif
