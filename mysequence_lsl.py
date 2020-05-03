###exercise1
# a
sequence="AAAACCTCTCTGTTCAGCACTTCCTCTCTCTTGGTCTGGTCTCAACGGTCACCATGGCGAGACCCTTGGAGGAGGCCCTGGATGTAATAGTGTCCACCTTCCACAAATACTCAGGCAACGAGGGTGACAAGTTCAAGCTGAACAAGACAGAGCTCAAGGAGCTACTGACCAGGGAGCTGCCTAGCTTCCTGGGGAGAAGGACAGACGAAGCTGCATTCCA"
dna=['A','C','G','T']
Dna=['Adenine','Cytosine','Guanine','Thymine']
print("Base counts:")
# for i in dna:
#     number_i=sequence.count(i,0,len(sequence))
#     print(i+":"+ str(number_i))
for i in Dna:
        number_i=sequence.count(i[0],0,len(sequence))
        print(i+":"+ str(number_i))

###exercise1
# b
# sequence="AAAACCTCTCTGTTCAGCACTTCCTCTCTCTTGGTCTGGTCTCAACGGTCACCATGGCGAGACCCTTGGAGGAGGCCCTGGATGTAATAGTGTCCACCTTCCACAAATACTCAGGCAACGAGGGTGACAAGTTCAAGCTGAACAAGACAGAGCTCAAGGAGCTACTGACCAGGGAGCTGCCTAGCTTCCTGGGGAGAAGGACAGACGAAGCTGCATTCCA"
your_string = input("Please write down your sequence and i will count the number of each base:")
Dna = ['Adenine', 'Cytosine', 'Guanine', 'Thymine']
def baseCounts(a):
    print("Base counts:")
    for i in Dna:
        number_i = your_string.count(i[0], 0, len(your_string))
        print(i + ":" + str(number_i))

baseCounts(your_string)


###exercise1
# c
# sequence="AAAACCTCTCTGTTCAGCACTTCCTCTCTCTTGGTCTGGTCTCAACGGTCACCATGGCGAGACCCTTGGAGGAGGCCCTGGATGTAATAGTGTCCACCTTCCACAAATACTCAGGCAACGAGGGTGACAAGTTCAAGCTGAACAAGACAGAGCTCAAGGAGCTACTGACCAGGGAGCTGCCTAGCTTCCTGGGGAGAAGGACAGACGAAGCTGCATTCCA"
your_string_format=input("Is your sequence DNA, RNA or amino acid?")
# your_string_format=your_string_format.lower
your_string=input("Please write down your sequence and i will count the number of each base:")
dna=list('ACGT')
rna=list('ACGU')
amino_acid=list('ACDEFGHIKLMNPQRSTVWY')
def baseCounts_dna(a):
    print("Base counts:")
    for i in dna:
        number_i=your_string.count(i,0,len(your_string))
        print(i+":"+ str(number_i))
def baseCounts_rna(a):
    print("Base counts:")
    for i in rna:
        number_i=your_string.count(i,0,len(your_string))
        print(i+":"+ str(number_i))
def baseCounts_amino_acid(a):
    print("Base counts:")
    for i in amino_acid:
        number_i=your_string.count(i,0,len(your_string))
        print(i+":"+ str(number_i))
if your_string_format=="DNA":
    print(baseCounts_dna(your_string))
elif your_string_format=="RNA":
    print(baseCounts_rna(your_string))
elif your_string_format=="amino acid":
    print(baseCounts_amino_acid(your_string))
###exercise1
# d
sequence = "AAAACCTCTCTGTTCAGCACTTCCTCTCTCTTGGTCTGGTCTCAACGGTCACCATGGCGAGACCCTTGGAGGAGGCCCTGGATGTAATAGTGTCCACCTTCCACAAATACTCAGGCAACGAGGGTGACAAGTTCAAGCTGAACAAGACAGAGCTCAAGGAGCTACTGACCAGGGAGCTGCCTAGCTTCCTGGGGAGAAGGACAGACGAAGCTGCATTCCA"
seDict={ }
Dna = ['Adenine', 'Cytosine', 'Guanine', 'Thymine']
def baseCounts(my_string):
    # print("Base counts:")
    for i in Dna:
        number_i = my_string.count(i[0], 0, len(my_string))
        seDict[i] = number_i

baseCounts(sequence)
print(seDict)


