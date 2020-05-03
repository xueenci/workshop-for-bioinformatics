###exercise 3
aDic = {}
Codon=[]
class Sequence():
    def __init__(self, sequence, alphabet, name='', gappy=False, annot=''):
        """ Construct a sequence from a string, an alphabet (gappy or not) and a name.
            The parameter gappy is for sequences when used in alignments. """
        for sym in sequence:
            if not sym in alphabet and (sym != '-' or not gappy):  # error check: bail out
                raise RuntimeError('Invalid symbol: ' + sym)
        self.sequence = sequence  # Store sequence
        self.alphabet = alphabet  # Store alphabet
        self.name = name  # Store name
        self.gappy = gappy
        self.annot = annot  # some annotation, e.g. species

    def __len__(self):  # the "len" operator
        return len(self.sequence)

    def __iter__(self):  # method that allows us to iterate over a sequence
        tsyms = tuple(self.sequence)
        return tsyms.__iter__()

    def __contains__(self, item):  # test for membership (the "in" operator)
        for sym in self.sequence:
            if sym == item:
                return True
        return False

    def __getitem__(self,
                    ndx):  # [ndx] operator (retrieve a specified index (or a "slice" of indices) of the sequence data.
        return self.sequence[ndx]

    def writeFasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        fasta = '>' + self.name + ' ' + self.annot + '\n'
        data = self.sequence
        nlines = (len(self.sequence) - 1) // 60 + 1
        for i in range(nlines):
            lineofseq = ''.join(data[i * 60: (i + 1) * 60]) + '\n'
            fasta += lineofseq
        return fasta

    def __str__(self):  # "pretty" print sequence
        str = self.name + ': '
        for sym in self:
            str += sym
        return str

    #     def count(self, findme):
    #         """ Get the number of occurrences of specified symbol """
    #         cnt = 0
    #         for sym in self.sequence:
    #             if findme == sym:
    #                 cnt = cnt + 1
    # #         return cnt
    #         aDic[findme]=cnt


    ###exercise 3a
    def count(self, findme):
        """Version modifyed: return a dictionary with the number of each base"""
        cnt = 0
        for sym in self.sequence:
            cnt = findme.count(sym, 0, len(findme))
            aDic[sym] = cnt


    ###exercise 3b
    se_2 = Sequence(
        "MHSSIVLATVLFVAIASASKTRELCMKSLEHAKVGTSKEAKQDGIDLYKHMFEHYPAMKKYFKHRENYTPADVQKDPFFIKQGQNILLACHVLCATYDDRETFDAYVGELMARHERDHVKVPNDVWNHFWEHFIEFLGSKTTLDEPTKHAWQEIGKEFSHEISHHGRHSVRDHCMNSLEYIAIGDKEHQKQNGIDLYKHMFEHYPHMRKAFKGRENFTKEDVQKDAFFVNKDTRFCWPFVCCDSSYDDEPTFDYFVDALMDRHIKDDIHLPQEQWHEFWKLFAEYLNEKSHQHLTEAEKHAWSTIGEDFAHEADKHAKAEKDHHEGEHKEEHH",
        Protein_Alphabet)
    se_2.count(
        "MHSSIVLATVLFVAIASASKTRELCMKSLEHAKVGTSKEAKQDGIDLYKHMFEHYPAMKKYFKHRENYTPADVQKDPFFIKQGQNILLACHVLCATYDDRETFDAYVGELMARHERDHVKVPNDVWNHFWEHFIEFLGSKTTLDEPTKHAWQEIGKEFSHEISHHGRHSVRDHCMNSLEYIAIGDKEHQKQNGIDLYKHMFEHYPHMRKAFKGRENFTKEDVQKDAFFVNKDTRFCWPFVCCDSSYDDEPTFDYFVDALMDRHIKDDIHLPQEQWHEFWKLFAEYLNEKSHQHLTEAEKHAWSTIGEDFAHEADKHAKAEKDHHEGEHKEEHH")
    print(aDic)


    ###exercise 3c
    def reverseComplement(self, sequence):
        """ reverse the sequence"""
        if self.alphabet == DNA_Alphabet:
            newSequence = ''
            for symbol in sequence[::-1]:
                #                 newSequence += symbol
                if symbol == "A":
                    newSequence += "T"
                elif symbol == "T":
                    newSequence += "A"
                elif symbol == "C":
                    newSequence += "G"
                elif symbol == "G":
                    newSequence += "C"
            return newSequence
        else:
            print("sorry, it is not a DNA sequence")

    ###exercise 3d
    def reverseComplement(self):
        """ reverse the sequence"""
        if self.alphabet == DNA_Alphabet:
            newSequence = ""
            for symbol in self.sequence[::-1]:
                #                 newSequence += symbol
                if symbol == "A":
                    newSequence += "T"
                elif symbol == "T":
                    newSequence += "A"
                elif symbol == "C":
                    newSequence += "G"
                elif symbol == "G":
                    newSequence += "C"
            return newSequence
        elif self.alphabet == RNA_Alphabet:
            newSequence = ""
            for symbol in self.sequence[::-1]:
                #                 newSequence += symbol
                if symbol == "A":
                    newSequence += "U"
                elif symbol == "T":
                    newSequence += "A"
                elif symbol == "C":
                    newSequence += "G"
                elif symbol == "G":
                    newSequence += "C"
        else:
            print("sorry, it is not a DNA or a RNA sequence")
    def find(self, findme):
        """ Find the position of the specified symbol or sub-sequence """
        return self.sequence.find(findme)

    # exercise 3e
    def translateDNA(self, rf, dr):
        if self.alphabet == DNA_Alphabet:
            i = rf
            while i < len(self.sequence):
                if dr == True:
                    Codon.append(self.sequence[i:i + 3])
                elif dr == False:
                    #                     self.sequence=self.reverseComplement()
                    newSequence = ""
                    for symbol in self.sequence[::-1]:
                        newSequence += symbol
                    Codon.append(newSequence[i:i + 3])
                i = i + 3
                for a in Codon:
                    if len(a) != 3:
                        Codon.remove(a)
            for i in Codon:
                for key, value in standardTranslation.items():
                    if i == key:
                        trans.append(value)
        else:
            print("sorry, it is not a DNA sequence")