#!/usr/bin/python
"""
Python module to map isolate mutations to likely drug resistances

Author: Ken Ramey
Date:   26 February 2018

CPTR Project - Critical Path Institute
"""
import os
import datetime
import logging
import argparse
from lxml import etree as ET
from sample_annotation import AnnotationList
from sample_coverage import SampleCoverage
from mutation_list import MutationList
from resistance_list import ResistanceList
from configparser import ConfigParser
from low_quals import LowQuals


class SurveillanceReport:
    """
    This class accumulates data for a surveillance report and stores it
    in an XML format for later processing.
    """
    def __init__(self, data_dir, sample_id):
        """
        Initialize the surveillance report data with local values and then
        get the required data from various pipeline output files.
        :param data_dir: Directory holding the input files for the report.
        :type data_dir: str
        :param sample_id: Identifier for the sample to be processed.
        :type sample_id: str
        """
        logging.debug('Initializing report object.')
        self.data_dir = data_dir
        self.sample_id = sample_id
        self.page = ET.Element('surveillance_report')
        self.doc = ET.ElementTree(self.page)
        title_elt = ET.SubElement(self.page, 'title')
        title_elt.text = 'Sample Surveillance Report'
        sample_id_elt = ET.SubElement(self.page, 'sample_id')
        sample_id_elt.text = sample_id
        date_elt = ET.SubElement(self.page, 'date')
        date_elt.text = str(datetime.datetime.now())
        pipeline_id_elt = ET.SubElement(self.page, 'pipeline')
        pipeline_name_elt = ET.SubElement(pipeline_id_elt, 'name')
        pipeline_name_elt.text = 'UVP'
        pipeline_version_elt = ET.SubElement(pipeline_id_elt, 'version')
        pipeline_version_elt.text = '1.1'
        self.lineage_elt = ET.SubElement(self.page, 'lineage')
        self.coverage_elt = ET.SubElement(self.page, 'coverage')
        self.coverage_gaps_elt = ET.SubElement(self.page, 'coverage_gaps')
        self.deletions_elt = ET.SubElement(self.page, 'deletions')
        self.mutations_elt = ET.SubElement(self.page, 'mutations')

    def report_lineage(self):
        """
        Get lineage data for isolate being processed.
        :return: Nothing
        :rtype: None
        """
        logging.debug('Determining isolate lineage.')
        f = open(os.path.join(self.data_dir, self.sample_id, self.sample_id + '.lineage_report.txt'), 'r')
        f.readline()
        fields = f.readline().strip().split('\t')
        lineage_code_elt = ET.SubElement(self.lineage_elt, 'code')
        lineage_code_elt.text = fields[3]
        lineage_text_elt = ET.SubElement(self.lineage_elt, 'name')
        lineage_text_elt.text = fields[2]
        f.close()

    def write_report_xml(self, xml_filename):
        """
        Write surveillance report data to file in xml format.
        :param xml_filename: Name of file to be created.
        :type xml_filename: str
        :return: Nothing
        :rtype: None
        """
        logging.debug('Creating XML version of data.')
        f = open(xml_filename, 'wb')
        self.doc.write(f)

    def add_coverage_section(self, coverage_map, genome_coverage):
        """
        Add details of region coverage to report.
        :param coverage_map: List of regions and their coverage amounts.
        :type coverage_map: dict
        :param genome_coverage: Coverage data for whole genome.
        :type genome_coverage: list
        :return: Nothing
        :rtype: None
        """
        logging.debug('Adding coverage data to report object.')
        gene_elt = ET.SubElement(self.coverage_elt, 'gene', attrib={'name': 'whole_genome'})
        gene_depth_elt = ET.SubElement(gene_elt, 'depth')
        gene_depth_elt.text = str(genome_coverage[0])
        gene_percent_elt = ET.SubElement(gene_elt, 'percent')
        gene_percent_elt.text = str(genome_coverage[1])
        for key, value in coverage_map.items():
            gene_elt = ET.SubElement(self.coverage_elt, 'gene', attrib={'name': key})
            gene_depth_elt = ET.SubElement(gene_elt, 'depth')
            gene_depth_elt.text = str(value[0])
            gene_percent_elt = ET.SubElement(gene_elt, 'percent')
            gene_percent_elt.text = str(value[1])

    def add_coverage_gaps_section(self, gap_map):
        """
        Document gaps in genome coverage by region.
        :param gap_map: List of regions with coverage details.
        :type gap_map: list
        :return: Nothing
        :rtype: None
        """
        logging.debug('Adding coverage gap data to report.')
        for gap in gap_map:
            gap_elt = ET.SubElement(self.coverage_gaps_elt, 'gene', attrib={'name': gap[0]})
            depth_elt = ET.SubElement(gap_elt, 'depth')
            depth_elt.text = str(gap[1][0])
            region_cov_elt = ET.SubElement(gap_elt, 'percent')
            region_cov_elt.text = str(gap[1][1])

    def add_low_quality_section(self, low_qual_list):
        """
        Create XML section that lists low quality sequence segments.
        :param low_qual_list: List of segments identified as low quality.
        :type low_qual_list: list
        :return: Nothing
        :rtype: None
        """
        low_qual_elt = ET.SubElement(self.page, 'low_quality')
        for segment in low_qual_list:
            segment_elt = ET.SubElement(low_qual_elt, 'segment', attrib={'refpos': segment[1]})
            ref_elt = ET.SubElement(segment_elt, 'ref')
            ref_elt.text = segment[2]
            alt_elt = ET.SubElement(segment_elt, 'alt')
            alt_elt.text = segment[3]
            qual_det_elt = ET.SubElement(segment_elt, 'qual_det')
            qual_det_elt.text = segment[4]

    def add_deletions_section(self, deletion_list):
        logging.debug('Adding deletions list to report.')
        for delete in deletion_list:
            deletion_elt = ET.SubElement(self.deletions_elt, 'loci', attrib={'name': delete[0]})
            del_type_elt = ET.SubElement(deletion_elt, 'type')
            del_type_elt.text = delete[1]

    def add_mutation_list_section(self, resistance_mutations):
        """
        Add resistance-significant mutations to report data.
        :param resistance_mutations: List of mutations along with drug data.
        :type resistance_mutations: list
        :return: Nothing
        :rtype: None
        """
        logging.debug('Adding drug resistance mutation data to report.')
        for resistance in resistance_mutations:
            gene_elt = ET.SubElement(self.mutations_elt, 'snp', attrib={'gene': resistance[0]})
            for index in range(1, len(resistance)):
                resistance_record = resistance[index]
                nuchange_elt = ET.SubElement(gene_elt, 'nuchange',
                                                attrib={'name': resistance_record[0], 'aachange': resistance_record[1]})
                annot_elt = ET.SubElement(nuchange_elt, 'annotation')
                annot_elt.text = resistance_record[5]
                cpos_elt = ET.SubElement(nuchange_elt, 'codonpos')
                cpos_elt.text = resistance_record[6]
                refpos_elt = ET.SubElement(nuchange_elt, 'refpos')
                refpos_elt.text = resistance_record[2]
                refnuc_elt = ET.SubElement(nuchange_elt, 'refnuc')
                refnuc_elt.text = resistance_record[3]
                altnuc_elt = ET.SubElement(nuchange_elt, 'altnuc')
                altnuc_elt.text = resistance_record[4]
                try:
                    for drug in resistance_record[7]:
                        resistance_elt = ET.SubElement(nuchange_elt, 'resistance')
                        drug_elt = ET.SubElement(resistance_elt, 'drug')
                        drug_elt.text = drug[0]
                        conf_elt = ET.SubElement(resistance_elt, 'confidence')
                        conf_elt.text = drug[1]
                except IndexError:
                    continue


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility to generate formatted report from UVP output.')
    parser.add_argument('-d', help='Directory holding UVP output.', required=True)
    parser.add_argument('-s', help='Name of sample to be processed.', required=True)
    parser.add_argument('-r', help='Reference list of mutations.', required=True)
    parser.add_argument('-o', help='Directory to hold XML output.', default='.')
    args = parser.parse_args()
    config = ConfigParser()
    config.read('reporting.cfg')
    defaults = config['Default']
    log_level = defaults['log_level'].strip()
    if log_level == 'DEBUG':
        lvl = logging.DEBUG
    elif log_level == 'INFO':
        lvl = logging.INFO
    else:
        lvl = logging.WARN
    logging.basicConfig(filename='surveillance.log', level=lvl)
    target_file = defaults['target_list'].strip()
    logging.info('Report creation at {}'.format(datetime.datetime.now()))
    data_directory = args.d
    reference_mutation_list = args.r
    output_directory = args.o
    sample_name = args.s
    # create reference mutation list
    reference_list = MutationList(reference_mutation_list)
    # read the pipeline annotation output
    pipeline_report = AnnotationList(data_directory, sample_name)
    report = SurveillanceReport(data_directory, sample_name)
    report.report_lineage()
    coverage = SampleCoverage(data_directory, sample_name, target_file)
    mutations = ResistanceList(data_directory, sample_name)
    lowquals = LowQuals(data_directory, sample_name)
    report.add_coverage_section(coverage.coverage_dict, [coverage.genome_depth, coverage.genome_percentage])
    report.add_coverage_gaps_section(coverage.coverage_gaps)
    report.add_deletions_section(coverage.deletions)
    report.add_mutation_list_section(mutations.resistance_list)
    report.add_low_quality_section(lowquals.low_quals)

    report.write_report_xml(sample_name + '_report.xml')
