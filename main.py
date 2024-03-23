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
    
    # Initialize NW class and specify gap pentalites 
    nw_62 = NeedlemanWunsch(sub_matrix_file= "./substitution_matrices/BLOSUM62.mat", gap_open= -10, gap_extend= -1)

    # perform alignments, add description and header, use tuples to make immutable
    gg_hs = (*nw_62.align(hs_seq, gg_seq), "Gallus gallus", gg_header)
    mm_hs = (*nw_62.align(hs_seq, mm_seq), "Mus musculus", mm_header)
    br_hs = (*nw_62.align(br_seq, hs_seq), "Balaeniceps rex", br_header)
    tt_hs = (*nw_62.align(tt_seq, hs_seq), "Tursinops trucatus", tt_header)

    # make list of alignments 
    aligns = [gg_hs, mm_hs, br_hs, tt_hs]

    # sort alignments based on scores
    sorted_aligns = sorted(aligns, key=lambda x: x[0], reverse=True)
    
    # generate list of aligned headers, print result
    headers_sorted = [sorted_aligns[i][3] for i in range(len(sorted_aligns))]
    print("Sorted in order of similarity to humans:")
    for i in headers_sorted:
        print("\t", i)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
        
    # generate list of sorted headers and scores, print result
    scores_sorted = [ (sorted_aligns[i][3], sorted_aligns[i][0]) for i in range(len(sorted_aligns))]
    print("\nAlignment scores between each species BRD2 and human BRD2:")
    for i in scores_sorted:
        print("\t", i)
    
    

if __name__ == "__main__":
    main()
