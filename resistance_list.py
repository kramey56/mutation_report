#!/usr/bin/python
"""
Python module to build the resistance list from the sample annotation and the graded mutation list

Author: Ken Ramey
Date:   28 February 2018

CPTR Project - Critical Path Institute
"""
from typing import List

from mutation_list import MutationList
from sample_annotation import AnnotationList


def confidence(lr):
    if lr >= 10.0:
        return 'High'
    elif 10 > lr >= 5:
        return 'Medium'
    else:
        return 'Low'


class ResistanceList:
    def __init__(self, data_dir, sample_id):
        self._resistance_list = []               # type: List
        reference_list = MutationList('Resistance_Rep_whocc_0218.csv')
        isolate_annotations = AnnotationList(data_dir, sample_id)
        self.references = reference_list.mutation_list
        self.targets = isolate_annotations.annotation_list
        gene_record = ['xxxx']                                          # type: List
        nuc_record = ['xxxx', 'xxxx', 'xxxx', 'xxxx', 'xxxx']           # type: List
        for target in self.targets:
            if target.gene != gene_record[0]:
                if len(nuc_record) > 7:
                    gene_record.append(nuc_record)
                    self._resistance_list.append(gene_record)
                gene_record = [target.gene]
                nuc_record = [target.nuchange, target.aachange, target.refpos, target.refnuc, target.altnuc,
                              target.annotation, target.codonpos]
            if target.nuchange != nuc_record[0]:
                gene_record.append(nuc_record)
                nuc_record = [target.nuchange, target.aachange, target.refpos, target.refnuc, target.altnuc,
                              target.annotation, target.codonpos]
            found, drug_list = self.find_drug_resistances(target.gene, target.nuchange)
            if found:
                nuc_record.append(drug_list)

    def find_drug_resistances(self, target_gene, target_nuchange):
        drug_resistances = []
        success = False
        for reference in self.references:
            if reference.gene == target_gene and reference.nuchange == target_nuchange:
                if len(drug_resistances) == 0:
                    drug_resistances = [[reference.drug, confidence(reference.likelihood)]]
                else:
                    drug_resistances.append([reference.drug, confidence(reference.likelihood)])
        if len(drug_resistances) > 0:
            success = True
        return success, drug_resistances

    @property
    def resistance_list(self):
        return self._resistance_list


if __name__ == '__main__':
    rl = ResistanceList('data', 'Rosenthal')
    for resistance in rl.resistance_list:
        print('length of {} record: {}'.format(resistance[0], len(resistance)))
        index = 0
        for item in resistance:
            print('item {} = {}'.format(index, item))
            index += 1
