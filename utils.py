import argparse


def parse_args():
    """
    Parse input arguments
    """
    parser = argparse.ArgumentParser(description='Binding Event Prediction')
    parser.add_argument('--ref', dest='ref', default='data/Homo_sapiens.GRCh38.dna.chromosome.1.fa', type=str)

    args = parser.parse_args()
    return args


def pwm():
    # position weight matrix
    A = [7, 15, 2, 1, 21, 0, 1, 0, 0, 1, 3]
    C = [1, 1, 9, 20, 0, 20, 0, 1, 0, 1, 5]
    G = [9, 4, 9, 0, 0, 0, 20, 0, 21, 18, 0]
    T = [4, 1, 1, 0, 0, 1, 0, 20, 0, 1, 13]
    len_fasta = len(A)
    return A, C, G, T, len_fasta


def PWM():
    A = [1038, 1019, 438, 44, 4600, 54, 169, 59, 23, 419, 961, 755]
    C = [1345, 1287, 3746, 4840, 59, 4525, 157, 162, 56, 2801, 1524, 1785]
    G = [1649, 1757, 472, 31, 170, 70, 4578, 63, 4817, 1198, 1036, 1348]
    T = [910, 879, 286, 27, 113, 293, 38, 4658, 46, 524, 1421, 1054]
    len_fasta = len(A)
    return A, C, G, T, len_fasta
