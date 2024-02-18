# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        # edit
        # print(dict_sub)
        # end edit
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing


        # mat_align_score = np.zeros_like()
        # mat_gaps = 
        # mat_backtrace =

        print(seqA)
        print(seqB)
        m = range(0, len(seqA)+1) # matrix rows
        n = range(0, len(seqB)+1) # matrix columns
        M = np.zeros( (len(seqA) + 1, len(seqB) + 1 ) )
        for i in m:
            for j in n:
                if i == 0:
                    M[i,j] = -np.inf
                if j == 0:
                    M[i,j] = -np.inf
                if i == 0 and j == 0:
                    M[i,j] = 0
        print(M)

        # highest gap penalty
        max_gap = self.gap_open * max(len(m), len(n))
        print("max gap:", max_gap)

        # extend penalty
        # self.gap_extend

        # gap penalty matrices
        I_m_y = np.zeros_like(M)
        for i in m:
            for j in n:
                if i == 0:
                    I_m_y[i,j] = self.gap_open * j # gap penalties along m columns (x-axis)
                if j == 0:
                    I_m_y[i,j] = -np.inf # infinite gap penalty 
                if i == 0 and j == 0:
                    I_m_y[i,j] = 0
        # print(I_m_y)

        I_n_x = np.zeros_like(M)
        for i in m:
            for j in n:
                if i == 0:
                    I_n_x[i,j] = -np.inf # infinite gap penalty 
                if j == 0:
                    I_n_x[i,j] = self.gap_open * i # gap penalties along n rows (y-axis)
                if j == 0:
                    if i == 0:
                        I_n_x[i,j] = 0
        # print("I_n_x: \n", I_n_x)
        # pass


        
        # TODO: Implement global alignment here

        # print(self.sub_dict.keys() )
        # print(self.sub_dict[('A', 'A')])
        # print(self.sub_dict.values())

        # consider each term in matrix M 
        for i in m:
            for j in n:
                if i - 1 in m and j - 1 in n: 
                    # current match-mismatch penalty
                    ij_score = self.sub_dict[(f'{seqA[i-1]}', f'{seqB[j-1]}')] 
                
                    # for the M matrix
                    M_i_1 = M[i - 1, j - 1] # get middle weights
                    # s = M[i,j]
                    M_score = M[i - 1, j - 1] + ij_score



                    # for the I_m_y (rows) matrix 
                    # if I_m_y[i , j - 1 ] < max_gap:
                    #     print("this is negative", (i,j) )
                    if I_m_y[i, j - 1] < max_gap:
                        # print("this is negative inf", (i,j) )
                        # I_m_y[i,j] = -np.inf
                        I_m_y_score = I_m_y[i , j - 1]
                    else:
                        I_m_y_score = I_m_y[i , j - 1] + self.gap_extend


                    # for the I_n_x (columns) matrix 
                    if I_n_x[i - 1, j] < max_gap:
                        # print("this is negative inf", (i,j) )
                        # I_n_x[i,j] = -np.inf
                        I_n_x_score = I_n_x[i - 1, j]
                    else:
                        I_n_x_score = I_n_x[i - 1, j] + self.gap_extend

                    # assign M[i,j] the highest score
                    M[i,j] = max( [M_score, I_m_y_score, I_n_x_score] )

                    I_n_x[i,j] = max( [M_score, I_n_x_score  ] )

                    I_m_y[i,j] = max( [M_score, I_m_y_score  ] )
                    print((i,j), M[i,j])
                # print(self.sub_dict[ (str(i),str(j))] )
        print("this is the M matrix \n",M)
        print(I_m_y)
        print(I_n_x)
        # pass      		
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # pass
        return (self.alignment_score, self.seqA_align, self.seqB_align), print("backtrace!")


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
