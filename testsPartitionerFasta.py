import os
import unittest
import fastaPartitionerIndex as fp
import random

import main

results = []
class TestPartitionOptions(unittest.TestCase):
    def setUp(self):
        self.path_data_file = './output_data/genes_index.fai'
        self.functions = fp.FunctionsFastaIndex(self.path_data_file)
        if not os.path.exists('./input_data/genes.fasta.fai'):
            main.generate_fasta_index_pyfaidx()
        with open('./input_data/genes.fasta.fai', "r") as f:
            self.data_assert = f.read().splitlines()

    def test_get_info_sequence(self):
        i = 0
        list = []
        for el in self.data_assert:
            if random.randint(0, 4) == 0:
                list.append(el.split('\t'))
                if i > 500:
                    break
                i += 1
        list.extend([['tr|IAmFalse'], ['']])

        for i, el in enumerate(list):
            info = self.functions.get_info_sequence(el[0])
            if info['length'] is not None and info['offset'] is not None and info['offset_head'] is not None:
                self.assertEqual([el[0], info['length'], info['offset']], [el[0], int(el[1]), int(el[2])])

    def test_get_range_sequence(self):
        results.append(f"Test 'test_get_range_sequence'")
        list = [[9643, 36109], [60411, 70827], [113494, 120950], [473060, 717220], [949179, 957690]]
        for i in list:
            sequences = self.functions.get_sequences_of_range(i[0], i[1])
            if sequences:
                seq = sequences[0].replace('\n','') + '\t\t-\t\t' + sequences[-1].replace('\n','')
                results.append(f"{i[0]}-{i[1]}: {seq}")
            else:
                results.append("Null")

    def test_index_generated(self):
        with open(self.functions.data_path, 'r') as index:
            sequence = index.readline()
            if sequence:
                for seq_assert in self.data_assert:
                    info_assert = seq_assert.split('\t')
                    info = sequence.split(' ')
                    length = int(info[3])
                    sequence = index.readline()
                    while sequence and info[0] == sequence.split(' ')[0]:
                        length += int(sequence.split(' ')[3])
                        sequence = index.readline()
                    self.assertEqual([info[0], str(length), info[2]], [info_assert[0], info_assert[1], info_assert[2]])
                    if not sequence:
                        break


