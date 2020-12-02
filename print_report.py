#!/usr/bin/python
"""
Python module to read xml report and produce printed output

Author: Ken Ramey
Date:   08 March 2018

CPTR Project - Critical Path Institute
"""
import os
import csv
import logging
import argparse
import datetime
from lxml import etree as ET
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER
from reportlab.lib.units import inch
from configparser import ConfigParser


class PrintableReport:
    """
    This class represents the printable version of the SNP report. It translates
    from XML to either text, PDF, or HTML.
    """
    def __init__(self, xml_string):
        """
        Initialize the data structures used for creating documents.
        :param xml_string: #The XML representation of the report data.
        :type xml_string:
        """
        self._root = ET.fromstring(xml_string)
        # self._root = self._tree.getroot()

        self._title = self._root.find('title').text
        self._test_date = self._root.find('date').text
        self._sample_id = self._root.find('sample_id').text
        self._pipeline_name = self._root.find('pipeline/name').text
        self._pipeline_version = self._root.find('pipeline/version').text
        self._isolate_lineage_code = self._root.find('lineage/code').text
        self._isolate_lineage_name = self._root.find('lineage/name').text
        self._coverage_map = []
        self.build_coverage_map()
        self._gap_map = []
        self._delete_list = []
        self.build_deletions_list()
        # self.build_gap_map()
        self._mutation_list = []
        self.build_mutation_list()

    def __repr__(self):
        """
        This function creates the plain text representation of the XML data structure.
        :return: A printable representation of the data.
        :rtype: str
        """
        format_str = '{:^80}\n\nRun Date: {}\nPipeline: {} {}\n\nSample ID: {}\nLineage: {}({})'
        header_str = format_str.format(self._title, self._test_date, self._pipeline_name, self._pipeline_version,
                                       self._sample_id, self._isolate_lineage_code, self._isolate_lineage_name)
        cov_fmt = '{:^60} {:>8} {:>8}\n'
        cov_hdr = cov_fmt.format('region', 'depth', 'percent')
        cov_str = ''
        for item in self._coverage_map:
            cov_str = cov_str + cov_fmt.format(item[0], item[1], item[2])
        cov_sec = '\n\nCoverage:\n' + cov_hdr + cov_str

        gap_str = ''
        for item in self._gap_map:
            gap_str = gap_str + cov_fmt.format(item[0], item[1], item[2])
        gap_sec = '\n\nCoverage Gaps:\n' + cov_hdr + gap_str

        mutation_fmt = '{:<15} {:<18} {:<18} {:<18} {:<10}\n'
        mutation_hdr = mutation_fmt.format('gene', 'nucleotide change', 'amino acid change',
                                           'drug resistance', 'confidence')
        mut_str = ''
        for item in self._mutation_list:
            mut_str = mut_str + mutation_fmt.format(item[0], item[1], item[2], item[3], item[4])
        mut_sec = '\n\nResistance List:\n' + mutation_hdr + mut_str

        return header_str + cov_sec + gap_sec + mut_sec

    def build_coverage_map(self):
        """
        Function to build a printable representation of the genome coverages by whole genome
        and by specific regions.
        :return: Nothing
        :rtype: None
        """
        cov_list = self._root.find('coverage')

        for region in cov_list:
            item_list = [region.attrib['name'], region.find('depth').text, region.find('percent').text]
            self._coverage_map.append(item_list)

    def build_gap_map(self):
        """
        This function builds a printable representation of the sections of the genome that are
        not covered adequately.
        :return: Nothing
        :rtype: None
        """
        gap_list = self._root.find('coverage_gaps')
        for gene in gap_list:
            item_list = [gene.attrib['name'], gene.find('depth').text, gene.find('percent').text]
            self._gap_map.append(item_list)

    def build_deletions_list(self):
        del_list = self._root.find('deletions')
        for loci in del_list:
            self._delete_list.append(loci.text)

    def build_mutation_list(self):
        """
        A function to build a representation of the SNPs believed to be related to drug resistance.
        :return: Nothing
        :rtype: None
        """
        for gene in self._root.iterfind('mutations/snp'):
            for nuchange in gene.findall('nuchange'):
                mutation = [gene.attrib['gene'], nuchange.attrib['name'], nuchange.attrib['aachange']]
                recorded = False
                for resistance in nuchange.findall('resistance'):
                    drug_resistance = [resistance.find('drug').text, resistance.find('confidence').text]
                    mutation += drug_resistance
                    self._mutation_list.append(mutation)
                    mutation = [gene.attrib['gene'], nuchange.attrib['name'], nuchange.attrib['aachange']]
                    recorded = True
                if not recorded:
                    self._mutation_list.append(mutation)

    @property
    def title(self):
        """
        Report title
        :return: Document title
        :rtype: str
        """
        return self._title

    @title.setter
    def title(self, new_title):
        self._title = new_title

    @property
    def test_date(self):
        """
        The date the UVP analysis was run.
        :return: Run date
        :rtype: str
        """
        return self._test_date

    @test_date.setter
    def test_date(self, new_date):
        self._test_date = new_date

    @property
    def sample_id(self):
        """
        ID of the sample that was analyzed.
        :return: Sample ID
        :rtype: str
        """
        return self._sample_id

    @sample_id.setter
    def sample_id(self, new_id):
        self._sample_id = new_id

    @property
    def pipeline_name(self):
        """
        Name of the pipeline used.
        :return: Pipeline name.
        :rtype: str
        """
        return self._pipeline_name

    @pipeline_name.setter
    def pipeline_name(self, new_name):
        self._pipeline_name = new_name

    @property
    def pipeline_version(self):
        """
        Version of the pipeline used.
        :return: Pipeline version.
        :rtype: None
        """
        return self._pipeline_version

    @pipeline_version.setter
    def pipeline_version(self, new_version):
        self._pipeline_version = new_version

    @property
    def isolate_lineage_name(self):
        """
        Name of the lineage identified.
        :return: Lineage name.
        :rtype: str
        """
        return self._isolate_lineage_name

    @isolate_lineage_name.setter
    def isolate_lineage_name(self, new_name):
        self._isolate_lineage_name = new_name

    @property
    def isolate_lineage_code(self):
        """
        Numeric code of the identified lineage.
        :return: Lineage code.
        :rtype: str
        """
        return self._isolate_lineage_code

    @isolate_lineage_code.setter
    def isolate_lineage_code(self, new_code):
        self._isolate_lineage_code = new_code

    @property
    def coverage_map(self):
        """
        The report coverage map.
        :return: Coverage map.
        :rtype: list
        """
        return self._coverage_map

    @coverage_map.setter
    def coverage_map(self, new_list):
        self._coverage_map.append(new_list)

    @property
    def gap_map(self):
        """
        The map of insufficient coverage.
        :return: Coverage gap map.
        :rtype: list
        """
        return self._gap_map

    @gap_map.setter
    def gap_map(self, new_list):
        self._gap_map.append(new_list)

    @property
    def mutation_list(self):
        """
        The list of resistance-related SNPs.
        :return: Mutation List.
        :rtype: list
        """
        return self._mutation_list

    @mutation_list.setter
    def mutation_list(self, new_item):
        self.mutation_list.append(new_item)

    def write_text_version(self, filename):
        """
        Write out a copy of the text version of the analysis report.
        :param filename: Name of file to hold report.
        :type filename: str
        :return: Nothing
        :rtype: None
        """
        f = open(filename, 'w')
        f.write(self.__repr__())
        f.close()

    @staticmethod
    def pdf_coord(x, y, height, unit=1):
        """
        A helper function to translate to PDF coordinate system.
        :param x: Horizontal position
        :type x: int
        :param y: Vertical position
        :type y: int
        :param height: Vertical size of page.
        :type height: int
        :param unit: size/unit for unit conversion. Defaults to 1.
        :type unit: int
        :return: Converted x, y pair.
        :rtype: tuple
        """
        x, y = x * unit, height - y * unit
        return x, y

    @staticmethod
    def page_number(current_canvas, document):
        """
        Automatically add page numbers to document.
        :param current_canvas: The drawing area being used.
        :type current_canvas: reportlab.pdfgen.canvas.Canvas
        :param document: The pdf document which will hold the text.
        :type document: SimpleDocTemplate
        :return: Noting
        :rtype: None
        """
        current_canvas.saveState()
        current_canvas.setFont('Helvetica', 10)
        current_canvas.drawString(3.75 * inch, 0.5 * inch, 'Page {}'.format(document.page))
        current_canvas.restoreState()

    def write_pdf_version(self, filename):
        """
        This is the main function that is called to create a PDF version of the report.
        :param filename: Name of the file to hold the PDF report.
        :type filename: str
        :return: Nothing
        :rtype: None
        """
        doc = SimpleDocTemplate(filename, pagesize=letter)
        styles = getSampleStyleSheet()
        styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
        styles.add(ParagraphStyle(name='Left', alignment=TA_LEFT))
        styles.add(ParagraphStyle(name='Right', alignmnet=TA_RIGHT))
        title = styles['Title']
        normal = styles['Normal']
        flowables = []

        ptext = self._title
        para = Paragraph(ptext, style=title)
        flowables.append(para)
        flowables.append(Spacer(1, 12))
        report_date = self._test_date.split('.')
        tbl_data = [[Paragraph('<b>Date:</b> {}'.format(report_date[0]), style=normal),
                     '',
                     Paragraph('<b>Pipeline:</b> {} {}'.format(self._pipeline_name,
                                                               self._pipeline_version), style=normal)],
                    [Paragraph('<b>Sample ID:</b> {}'.format(self._sample_id), style=normal),
                     '',
                     Paragraph('<b>Lineage:</b> {}({})'.format(self._isolate_lineage_code, self._isolate_lineage_name),
                               style=normal)]
                    ]

        tbl = Table(tbl_data, 2.15 * inch)
        flowables.append(tbl)

        ptext = 'Coverage:'
        para = Paragraph(ptext, style=styles['Heading3'])
        flowables.append(para)

        cov_data = [['gene', 'depth', 'coverage %']]
        for row in self._coverage_map:
            cov_data.append(row)
        tbl_style = TableStyle([('ALIGNMENT', (-1, 0), (-1, -1), 'RIGHT'),
                                ('FONT', (0, 0), (-1, 0), 'Helvetica-Bold')])

        cov_tbl = Table(cov_data, 1.5 * inch)
        cov_tbl.setStyle(tbl_style)
        flowables.append(cov_tbl)
        flowables.append(Spacer(1, 12))

        """
        ptext = 'Coverage Gaps:'
        para = Paragraph(ptext, style=styles['Heading3'])
        flowables.append(para)
        gap_data = [['gene', 'depth', 'coverage %']]
        for row in self._gap_map:
            gap_data.append(row)
        gap_tbl = Table(gap_data, colWidths=[4.25 * inch, inch, inch], repeatRows=1)
        gap_tbl.setStyle(tblstyle=tbl_style)
        flowables.append(gap_tbl)
        flowables.append(Spacer(1, 12))
        """

        ptext = 'Deletions:'
        para = Paragraph(ptext, style=styles['Heading3'])
        flowables.append(para)
        if len(self._delete_list) > 0:
            del_data = [['gene', 'type']]
            for loci in self._delete_list:
                del_data.append(loci)
            del_tbl = Table(del_data, colWidths=(4.25 * inch))
            del_tbl.setStyle(tbl_style)
            flowables.append(del_tbl)
        else:
            ptext = '     None'
            para = Paragraph(ptext, style=ParagraphStyle('default'))
            flowables.append(para)
        flowables.append(Spacer(1, 12))

        ptext = 'Resistance List'
        para = Paragraph(ptext, style=styles['Heading3'])
        flowables.append(para)
        mut_data = [['gene', 'nucleotide change', 'amino acid change', 'drug', 'confidence']]
        for gene in self._mutation_list:
            mut_data.append(gene)
        mut_tbl = Table(mut_data, repeatRows=1)
        mut_tbl.setStyle(tbl_style)
        flowables.append(mut_tbl)

        doc.build(flowables, onLaterPages=self.page_number)

    def write_csv_data(self, csv_filename):
        fields = ['country_code', 'drs_code', 'nrl_code', 'srl_code', 'culture_nrl', 'culture_srl',
                  'cluster', 'smear', 'xpert_tb', 'gyra_nt_sanger', 'gyra_aa_sanger', 'gyra_nt_ngs',
                  'gyra_aa_ngs', 'gyra_nt_wgs', 'gyra_aa_wgs', 'gyra_mut', 'gyra_mut_name',
                  'gyra_mut_class_ofx_lfx', 'gyra_mut_class_mfx', 'gyrb_nt_sanger', 'gyrb_aa_sanger',
                  'gyrb_nt_ngs', 'gyrb_aa_ngs', 'gyrb_nt_wgs', 'gyrb_aa_wgs', 'gyrb_mut', 'gyrb_mut_name',
                  'gyrb_mut_class_ofx_lfx', 'gyrb_mut_class_mfx', 'pnca_nt_sanger', 'pnca_aa_sanger',
                  'pnca_nt_ngs', 'pnca_aa_ngs', 'pnca_nt_wgs', 'pnca_aa_wgs', 'pnca_mut', 'pnca_mut_name',
                  'pnca_mut_class', 'rpob_nt_sanger', 'rpob_aa_sanger', 'rpob_nt_ngs', 'rpob_aa_ngs',
                  'rpob_nt_wgs', 'rpob_aa_wgs', 'rpob_mut', 'rpob_mut_name', 'rpob_mut_class',
                  'inha_nt_sanger', 'inha_aa_sanger', 'inha_nt_ngs', 'inha_aa_ngs', 'inha_nt_wgs',
                  'inha_aa_wgs', 'inha_mut', 'inha_mut_name', 'inha_mut_class', 'katg_nt_sanger', 'katg_aa_sanger',
                  'katg_nt_ngs', 'katg_aa_ngs', 'katg_mut', 'katg_mut_name', 'katg_mut_class', 'rrs_nt_sanger',
                  'rrs_aa_sanger', 'rrs_nt_ngs', 'rrs_aa_ngs', 'rrs_nt_wgs', 'rrs_aa_wgs', 'rrs_mut_sli',
                  'rrs_mut_name_sli', 'rrs_mut_sm', 'rrs_mut_class_kan', 'rrs_mut_class_amk', 'dr', 'eis_nt_sanger',
                  'eis_aa_sanger', 'eis_nt_ngs', 'eis_aa_ngs', 'eis_nt_wgs', 'eis_aa_wgs', 'eis_mut', 'eis_mut_name',
                  'eis_mut_class_kan', 'rif', 'rif_pheno', 'rif_pheno_type', 'rif_molecular', 'inh', 'inh_type',
                  'pza', 'waynes_test', 'ofx2_1', 'ofx2_2', 'ofx_final', 'mfx05', 'mfx2', 'lfx15', 'gfx2', 'sm',
                  'emb', 'cap', 'kan', 'amk', 'sex', 'age', 'agegrp', 'hiv', 'previous_treat', 'cat_previous_treat',
                  'treat_regimen', 'treat_duration', 'treat_outcoe', 'lineage', 'genome_mean_coverage',
                  'pnca_gene_meancoverage', 'gyra_gene_meancoverage', 'gyrb_gene_meancoverage', 'rpob_mean_coverage',
                  'katg_mean_coverage', 'inha_mean_coverage', 'rrs_mean_coverage', 'eis_mean_coverage', 'notes',
                  'excluded']
        who_goi = ['gyra', 'gyrb', 'pnca', 'rpob', 'inha', 'katg', 'rrs', 'eis']
        report_row = []
        for index in range(128):
            report_row.append('')
        csv_out = open(csv_filename, mode='w')
        writer = csv.writer(csv_out)
        writer.writerow(fields)
        for gene in self._root.iterfind('mutations/gene'):
            if gene.attrib['name'] in who_goi:
                for nuchange in gene.findall('nuchange'):
                    for resistance in nuchange.findall('resistance'):
                        drug_resistance = [gene.attrib['name'], nuchange.attrib['name'], nuchange.attrib['aachange'],
                                           resistance.find('drug').text, resistance.find('confidence').text]
                        self.mutation_list = drug_resistance


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility to generate formatted report from UVP output.')
    parser.add_argument('data_file', help='XML file to be processed.')
    parser.add_argument('-o', help='Directory to hold report output.', default='.')
    parser.add_argument('-f', help='Format of generated report.', choices=['HTML', 'PDF', 'TEXT'], default='HTML')
    parser.add_argument('-s', help='Sample name for report naming.', required=True)
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
    logging.info('Report creation at {}'.format(datetime.datetime.now()))
    output_directory = args.o
    sample_name = args.s
    data_file = args.data_file
    with open(data_file) as xml_file:
        xml_content = xml_file.read().replace('\n', ' ')
    surv_report = PrintableReport(xml_content)

    if args.f == 'HTML':
        dom = ET.fromstring(xml_content)
        xslt = ET.parse("html_report.xsl")
        transform = ET.XSLT(xslt)
        newdom = transform(dom)
        report_file = os.path.join(output_directory, sample_name + '_report.html')
        newdom.write(report_file, pretty_print=True)
    elif args.f == 'PDF':
        report_file = os.path.join(output_directory, sample_name + '_report.pdf')
        surv_report = PrintableReport(xml_content)
        surv_report.write_pdf_version(report_file)
    elif args.f == 'TEXT':
        report_file = os.path.join(output_directory, sample_name + '_report.txt')
        surv_report.write_text_version(report_file)
