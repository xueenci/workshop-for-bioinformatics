from guide import *

seqs = readFastaFile("sarsCov2.fasta", DNA_Alphabet)
seq = seqs[0]
# print(seq)
seq=Sequence(seq,DNA_Alphabet)
# print(seq)
# seq.reverseComplement()
# # # print(seq)
a=seq.reverseComplement()
a[501:520]
