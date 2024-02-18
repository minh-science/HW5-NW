# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix


    # TESTING    
    t1_seq, t1_header = read_fasta("./data/test_seq1.fa")
    t2_seq, t2_header = read_fasta("./data/test_seq2.fa")
    t3_seq, t3_header = read_fasta("./data/test_seq3.fa")
    t4_seq, t4_header = read_fasta("./data/test_seq4.fa")

    B = "WARAW"
    A = "WAAAW"
    
    nw_62 = NeedlemanWunsch(sub_matrix_file= "./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    # nw_62.align(A, B)
    # nw_62.align(t2_seq, t1_seq)
    # nw_62.align(t1_seq, t2_seq)
    nw_62.align(t4_seq, t3_seq)
    # end TESTING

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    # pass
    
    

if __name__ == "__main__":
    main()
