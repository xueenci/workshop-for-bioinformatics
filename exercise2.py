from guide import readFastaFile
seqs = readFastaFile("sarsCov2.fasta", DNA_Alphabet)
seq = seqs[0]
seq=Sequence(seq,DNA_Alphabet)
# print(seq)
seq.reverseComplement()