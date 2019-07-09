from utils import parse_args, PWM


class MYC:
    """
    An algorithm that parse a fasta file, calculates the consensus sequence of the motif
    from the PWM and finds the most likely position in the reference genome that binds to.
    """

    def __init__(self, reference, A, C, G, T, len_fasta):
        self.consensus_pattern = ''
        self.find_consensus(A, C, G, T, len_fasta)
        ref = self.parse_fasta(reference)
        self.prediction = self.longest_common_substring(ref, self.consensus_pattern)

    def longest_common_substring(self, reference, pattern):
        """
        Dynamic programming algorithms to find the most likely binding site in the human genome based based on the lenght og the occurence
        :param reference: The human genome chromosome 1
        :param pattern: the trascription factor motif to be predicted
        :return: The most likely sequence in the first 100k bases to be bound by MYC
        """
        find_match = ""
        max_length = 0
        i_max, j_max = 0, 0
        w, h = len(pattern), len(reference)
        matrix = [[0 for x in range(w)] for y in range(h)]
        for i, row in enumerate(matrix):
            for j, item in enumerate(row):
                match = reference[i] == pattern[j]
                if match:
                    if i == 0 or j == 0:
                        matrix[i][j] = 1
                    else:
                        string_len = matrix[i - 1][j - 1] + 1
                        if string_len > max_length:
                            max_length = string_len
                            i_max, j_max = i, j
                        matrix[i][j] = string_len
        if max_length > 0:
            for i in range(1, max_length + 1):
                if i == 1:
                    print(
                        f'The most likely binding site is at position <{i_max - max_length + i +1}> in the reference genome')
                find_match = find_match + reference[i_max - max_length + i]

        return find_match

    def find_consensus(self, A, C, G, T, len_fasta):
        output_num = [0 for i in range(len_fasta)]

        for cons_max in range(len_fasta):
            output_num[cons_max] = max(A[cons_max], C[cons_max], G[cons_max], T[cons_max])
            if A[cons_max] >= C[cons_max] and A[cons_max] >= G[cons_max] and A[cons_max] >= T[cons_max]:
                self.consensus_pattern += 'A'

            elif C[cons_max] >= A[cons_max] and C[cons_max] >= G[cons_max] and C[cons_max] >= T[cons_max]:
                self.consensus_pattern += 'C'

            elif G[cons_max] >= A[cons_max] and G[cons_max] >= C[cons_max] and G[cons_max] >= T[cons_max]:
                self.consensus_pattern += 'G'

            elif T[cons_max] >= C[cons_max] and T[cons_max] >= G[cons_max] and T[cons_max] >= A[cons_max]:
                self.consensus_pattern += 'T'

        assert len(self.consensus_pattern) == len_fasta
        # print(len(output_string), len_fasta)

    def parse_fasta(self, ref):
        fin = open(ref, 'rb')
        from itertools import groupby
        faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))
        for header in faiter:
            header_fasta = str(header.__next__(), 'utf-8')
            long_name = header_fasta.strip().replace('>', '')
            name = long_name.split()[0]
            sequence = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
            # print(name)
            # print(f'what is found normal: {seq[35075:35084]}')
            # print(f'what is found reverse: {seq[35076:35086]}')
            return sequence[:100000]


if __name__ == '__main__':
    args = parse_args()
    print('Called with args:')
    print(args)
    reference = args.ref
    print('##############################################################################')
    out = MYC(reference, PWM()[0], PWM()[1], PWM()[2], PWM()[3], PWM()[4])
    print(f'The predicted binding site in the ref genome: {out.prediction}')
    print(f'The original pattern to be found: {out.consensus_pattern}')
    print('##############################################################################')
    print(out.consensus_pattern)
    print('  ||||||||')
    print(f'  {out.prediction}')
