#!/usr/bin/python
"""
Python module build the reference mutations list

Author: Ken Ramey
Date:   26 February 2018

CPTR Project - Critical Path Institute
"""
from collections import namedtuple


class MutationList:
    # Get reference list of graded mutations and map infinite values to 9999
    def __init__(self, reference_file):
        Mutation = namedtuple('Mutation', ['gene', 'nuchange', 'aachange', 'drug', 'pvalue', 'likelihood'])
        self.mutations = []
        self.item_count = 0
        with open(reference_file) as f:
            f.readline()
            f.readline()
            for line in f:
                fields = line.strip().split(',')
                graded_mutation = Mutation(fields[1], fields[6], fields[7], fields[0],
                                           float(fields[22]), float(fields[17]))
                if graded_mutation.pvalue < 0.05:
                    self.mutations.append(graded_mutation)
                    self.item_count += 1

    @property
    def size(self):
        return self.item_count

    @property
    def mutation_list(self):
        return sorted(self.mutations)


if __name__ == '__main__':
    mlist = MutationList('Resistance_Rep_whocc_0218.csv')
    print('Items in mutations list: {}'.format(mlist.size))
    for m in mlist.mutation_list:
        print(m)
