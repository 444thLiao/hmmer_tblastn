from Bio import SeqIO
aa_dict = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S",
           "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
           "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
           "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
           "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
           "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
           "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N",
           "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
           "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
           "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E",
           "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

def translate(offset, seq):
    """Does the actual translating"""
    tmp_str = ''
    for x in range(int(offset), len(seq), 3):
        codon = seq[x:x + 3]
        if len(codon) < 3:
            break
        elif "N" in codon:
            tmp_str += "X"
        else:
            tmp_str += "%s" % aa_dict.get(codon,'X')
    return tmp_str

def generate_six_frames(fasta):
    fasta = SeqIO.parse(fasta,format='fasta')
    fasta = list(fasta)
    seqs = [(_.description,str(_.seq).upper()) for _ in fasta]
    seq_rc = [(_.description,str(_.reverse_complement().seq).upper()) for _ in fasta]

    each_frame = [(_,i,translate(i,v)) for _,v in seqs for i in range(3)] + \
                 [(_, 'rev'+str(i), translate(i, v)) for _, v in seq_rc for i in range(3)]
    return each_frame


