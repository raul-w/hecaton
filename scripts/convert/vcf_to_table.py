#!/usr/bin/env python3

"""
Convert VCF files to tabular file, script based on VariantsToTable of GATK
"""

import argparse
import pysam
import logging
import sys


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Convert VCF file to tabular file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--vcf_fn", type=str,
                        help="File containing path to VCF file")
    parser.add_argument("-f", "--fields", type=str, nargs='*',
                        help="Standard VCF or INFO fields that will be included in the tabular file")
    parser.add_argument("-gf", "--genotype_fields", type=str, nargs='*',
                        help="Genotype fields that will be included in the tabular file")
    parser.add_argument("-o", "--output_fn", type=str,
                        help="Name of tabular output file")
    args = parser.parse_args(in_args)
    return args

def vcf_to_table(input_fn, output_fn, fields=None, genotype_fields=None):
    """
    Write specific fields of a VCF file to a tabular format

    :param input_fn: Path to VCF file
    :param output_fn: Name of tabular output file
    :param fields: Standard VCF or INFO fields that will be included in the tabular file
    :param genotype_fields: Genotype fields that will be included in the tabular file
    :return: 0 (integer)
    """
    if not fields and not genotype_fields:
        raise ValueError("At least one standard, INFO, or genotype field needs to be provided as input")
    with pysam.VariantFile(input_fn) as vcf, open(output_fn, 'w') as output_file:
        # create header
        samples = list(vcf.header.samples)
        samples.sort()
        header_fields = []
        standard_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "END", "HOM-VAR", "VAR", "SAMPLES-VAR", "DEL_SUPPORTED"]
        standard_fields.extend(vcf.header.info)
        if fields:
            for field in fields:
                if field not in standard_fields:
                    raise ValueError(
                        "{} not defined as standard or INFO field".format(field))
                if field == "DEL_SUPPORTED":
                    header_fields.extend(["DEL_SUPPORTED", "DEL_UNSUPPORTED", "NON_DEL_SUPPORTED", "NON_DEL_UNSUPPORTED"])
                else:
                    header_fields.append(field)
        if genotype_fields:
            for genotype_field in genotype_fields:
                if genotype_field not in vcf.header.formats:
                    raise ValueError(
                        "{} not defined as genotype field".format(genotype_field))
                else:
                    for sample in samples:
                        header_fields.append(sample + "." + genotype_field)
        # write header to output file
        header = '\t'.join(header_fields)
        output_file.write(header)
        output_file.write("\n")
        # extract all fields
        for record in vcf.fetch():
            record_fields = []
            if fields:
                for field in fields:
                    if field == "CHROM":
                        record_fields.append(str(record.chrom))
                    elif field == "POS":
                        record_fields.append(str(record.pos))
                    elif field == "ID":
                        record_fields.append(str(record.id))
                    elif field == "REF":
                        record_fields.append(str(record.ref))
                    elif field == "ALT":
                        alts = []
                        for alt in record.alts:
                            alts.append(str(alt))
                        alts = ",".join(alts)
                        record_fields.append(alts)
                    elif field == "QUAL":
                        record_fields.append(str(record.qual))
                    elif field == "FILTER":
                        record_fields.append(str(record.filter.keys()[0]))
                    elif field == "END":
                        record_fields.append(str(record.stop))
                    elif field == "HOM-VAR":
                        hom_calls = 0
                        # compute total number of homozygous variant calls
                        if "GT" not in record.format:
                            logging.warning("GT not in format, HOM-VAR will be NA")
                            record_fields.append("NA")
                        else:
                            for sample in samples:
                                genotype = record.samples[sample]["GT"]
                                if genotype == (1, 1):
                                    hom_calls += 1
                            record_fields.append(str(hom_calls))
                    elif field == "VAR":
                        var_calls = 0
                        # compute total number of variant calls
                        if "GT" not in record.format:
                            logging.warning("GT not in format, HOM-VAR will be NA")
                            record_fields.append("NA")
                        else:
                            for sample in samples:
                                genotype = record.samples[sample]["GT"]
                                non_variants = [(0, 0), (None, None), (None, 0)]
                                if genotype not in non_variants:
                                    var_calls += 1
                            record_fields.append(str(var_calls))
                    elif field == "SAMPLES-VAR":
                        # get list of samples that contain a non-ref variant
                        samples_var = []
                        # compute total number of variant calls
                        if "GT" not in record.format:
                            logging.warning("GT not in format, SAMPLES-VAR will be NA")
                            record_fields.append("NA")
                        else:
                            for sample in samples:
                                genotype = record.samples[sample]["GT"]
                                non_variants = [(0, 0), (None, None), (None, 0)]
                                if genotype not in non_variants:
                                    samples_var.append((sample, genotype))
                            record_fields.append(str(samples_var))
                    elif field == "DEL_SUPPORTED":
                        del_supported_calls = 0
                        del_unsupported_calls = 0
                        nondel_supported_calls = 0
                        nondel_unsupported_calls = 0
                        if record.info["SVTYPE"] != "DEL":
                            record_fields.extend(["NA", "NA", "NA", "NA"])
                        else:
                            for sample in samples:
                                genotype = record.samples[sample]["GT"]
                                dhffc = record.samples[sample]["DHFFC"]
                                non_variants = [(0, 0), (None, None),
                                                (None, 0)]
                                if genotype in non_variants:
                                    if dhffc >= 0.75:
                                        nondel_supported_calls += 1
                                    else:
                                        nondel_unsupported_calls += 1
                                elif genotype == (0, 1):
                                    if dhffc >= 0.25 and dhffc < 0.75:
                                        del_supported_calls += 1
                                    else:
                                        del_unsupported_calls += 1
                                elif genotype == (1, 1):
                                    if dhffc >= 0 and dhffc < 0.25:
                                        del_supported_calls += 1
                                    else:
                                        del_unsupported_calls += 1
                            record_fields.append(str(del_supported_calls))
                            record_fields.append(str(del_unsupported_calls))
                            record_fields.append(str(nondel_supported_calls))
                            record_fields.append(str(nondel_unsupported_calls))
                    elif field in record.info:
                        record_fields.append(str(record.info[field]))
                    else:
                        record_fields.append("NA")
            if genotype_fields:
                for genotype_field in genotype_fields:
                    for sample in samples:
                        format_field = record.samples[sample][genotype_field]
                        if type(format_field) == float:
                            format_field = round(format_field, 2)
                        record_fields.append(str(format_field))
            # write record field to tabular output
            record_line = "\t".join(record_fields)
            output_file.write(record_line)
            output_file.write("\n")
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    if not args.fields and not args.genotype_fields:
        raise ValueError("At least one standard, INFO, or genotype field needs to be provided as input")
    # write VCF fields to tabular output
    vcf_to_table(args.vcf_fn, args.output_fn, fields=args.fields,
                 genotype_fields=args.genotype_fields)

if __name__ == "__main__":
    main()
