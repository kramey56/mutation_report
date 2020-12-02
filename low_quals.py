#!/usr/bin/python
"""
Python module to read low_quals text file and extract data
for specific sample.

Author: Ken Ramey
Date:   07 May 2018

CPTR Project - Critical Path Institute
"""

import os


class LowQuals:
    def __init__(self, location, sample):
        """
        Create a list of low quality sequence segments.
        :param location: Directory holding low quality report.
        :type location: str
        :param sample: Name of low quality report file.
        :type sample: str
        """
        self._lowqual_list = []
        lowqual_file = os.path.join(location, sample, 'low_quals.txt')
        with open(lowqual_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if fields[0] == sample:
                    lowqual_rec = [fields[0], fields[1], fields[2], fields[3], fields[4]]
                    self._lowqual_list.append(lowqual_rec)

    @property
    def low_quals(self):
        return self._lowqual_list


if __name__ == '__main__':
    sample_name = 'Rosenthal_S81'
    low_qual_dir = 'C:\Users\kramey\PycharmProjects\mutation_report'
    lq = LowQuals(low_qual_dir, sample_name)
    for record in lq.low_quals:
        print(record)
