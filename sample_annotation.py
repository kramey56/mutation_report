#!/usr/bin/python
"""
Python module to read the annotation output for an isolate

Author: Ken Ramey
Date:   26 February 2018

CPTR Project - Critical Path Institute
"""
import os
from collections import namedtuple


class AnnotationList:
    # Get annotated results from pipeline
    def __init__(self, data_dir, sample_id):
        results_file = os.path.join(data_dir, sample_id, sample_id + '_Resistance_Final_annotation.txt')
        IsolateMutation = namedtuple('IsolateMutation', ['gene', 'nuchange', 'aachange', 'refpos', 'refnuc',
                                                         'altnuc', 'annotation', 'codonpos'])
        self.mutations = []
        self.item_count = 0

        with open(results_file, 'r') as f:
            f.readline()
            for line in f:
                fields = line.strip().split('\t')
                isolated_mutation = IsolateMutation(fields[16], fields[10][2:], fields[12][2:], fields[2], fields[3],
                                                    fields[4], fields[8], fields[15])
                self.mutations.append(isolated_mutation)
                self.item_count += 1

    @property
    def annotation_list(self):
        return sorted(self.mutations)

    @property
    def size(self):
        return self.item_count


if __name__ == '__main__':
    test_list = AnnotationList('data', 'Rosenthal')
    print('Size of list = {}'.format(test_list.size))
    print('List:')
    for m in test_list.annotation_list:
        print(m)
