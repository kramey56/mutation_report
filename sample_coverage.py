#!/usr/bin/python
"""
Python module to read the sample coverage output for an isolate

Author: Ken Ramey
Date:   28 February 2018

CPTR Project - Critical Path Institute
"""
import os


class SampleCoverage:
    def __init__(self, data_dir, sample_id, loci_of_interest):
        self.sample_coverage_file = os.path.join(data_dir, sample_id, sample_id + '_Coverage.txt')
        self.region_coverage_file = os.path.join(data_dir, sample_id, sample_id + '_genome_region_coverage.txt')
        self.deletions_file = os.path.join(data_dir, sample_id, sample_id + '_deleted_loci.txt')
        self._genome_depth = 0
        self.genome_percentage = 0
        with open(loci_of_interest, 'r') as f:
            line = f.read()
        self._genes_of_interest = line.strip().split(',')
        self._coverage_dict = {}
        self._coverage_gaps = []
        self._deletions = []
        self.get_sample_stats()
        self.get_region_stats()
        self.get_deletions()

    def get_sample_stats(self):
        f = open(self.sample_coverage_file, 'r')
        line = f.readline().strip()
        fields = line.split(':')
        self._genome_depth = fields[1]
        line = f.readline().strip()
        fields = line.split(':')
        self.genome_percentage = fields[1]
        f.close()

    def get_region_stats(self):
        with open(self.region_coverage_file, 'r') as f:
            f.readline()
            for line in f:
                fields = line.split('\t')
                if fields[3] in self._genes_of_interest:
                    self._coverage_dict[fields[3]] = [float(fields[5]), float(fields[6])]
                if float(fields[6]) < 90.0:
                    gap_record = [fields[3], [float(fields[5]), float(fields[6])]]
                    self._coverage_gaps.append(gap_record)

    def get_deletions(self):
        with open(self.deletions_file, 'r') as f:
            f.readline()
            for line in f:
                fields = line.split('\t')
                if fields[15] in self._genes_of_interest:
                    del_record = [fields[15], fields[8]]
                    self._deletions.append(del_record)

    @property
    def coverage_dict(self):
        return self._coverage_dict

    @property
    def coverage_gaps(self):
        return self._coverage_gaps

    @property
    def deletions(self):
        return self._deletions

    @property
    def genome_depth(self):
        return self._genome_depth


if __name__ == '__main__':
    t = SampleCoverage('data', 'Rosenthal', 'Combined_Targets.csv')
    print('Done')
