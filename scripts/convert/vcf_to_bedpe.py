#!/usr/bin/env python3

"""
Convert VCF files to BEDPE
"""

import argparse
import pysam
import sys

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Parse VCF files, extracting simple CNVs in BEDPE format"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_vcf", type=str,
                        help="Path to VCF file")
    parser.add_argument("-o", "--output_bedpe", type=str,
                        help="Name of BEDPE output file")
    parser.add_argument("-t", "--tool", type=str,
                        help="Tool used to call variants")
    args = parser.parse_args(in_args)
    return args

def process_vcf(input_vcf_fn, output_bedpe_fn, tool):
    """Process vcfs, extracting simple cnvs

    :param input_vcf_fn: Filename of input VCF file
    :param output_bedpe_fn: Filename of output BEDPE file
    :param tool: Name of tool used to call variants
    :return: 0 (integer)
    """
    # collect DELS and DUPs vcf file and write to BEDPE
    with pysam.VariantFile(input_vcf_fn) as vcf, open(output_bedpe_fn, "w") as bedpe:
        # write header to bedpe
        header = ["#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B",
                  "END_B", "ID", "QUAL", "STRAND_A", "STRAND_B", "TYPE", "INSERTED_SEQUENCE", "TOOL",
                  "READ_PAIRS_DEFINED", "READ_PAIRS", "SPLIT_READS_DEFINED", "SPLIT_READS", "SIZE"]
        header = '\t'.join(header) + '\n'
        bedpe.write(header)
        # loop through each line of the vcf and write it to bedpe
        allowed_types = ["DEL", "DUP", "INS"]
        for record in vcf.fetch():
            if record.info["SVTYPE"] not in allowed_types:
                continue
            inserted_sequence = "."
            if record.info["SVTYPE"] == "INS":
                # get inserted sequence, if present
                inserted_sequence = "N"
                # get inserted sequence, if present
                alts = record.alts
                alt_field = alts[0]
                if alt_field != "<INS>":
                    inserted_sequence = alt_field[1:]
            # get read pair and split read information, if present
            read_pairs_defined = 0
            read_pairs = 0
            split_reads_defined = 0
            split_reads = 0
            # LUMPY format
            if "PE" in record.format:
                read_pairs_defined = 1
                read_pairs = record.samples[0]["PE"]
            # Manta format
            elif "PR" in record.format:
                read_pairs_defined = 1
                read_pairs = record.samples[0]["PR"][1]
            # DELLY format
            elif "PE" in record.info:
                read_pairs_defined = 1
                read_pairs = record.info["PE"]
            # LUMPY and Manta format
            if "SR" in record.format:
                # Manta format
                if isinstance(record.samples[0]["SR"], tuple):
                    split_reads_defined = 1
                    split_reads = record.samples[0]["SR"][1]
                # LUMPY format
                else:
                    split_reads_defined = 1
                    split_reads = record.samples[0]["SR"]
            # Delly format
            elif "SR" in record.info:
                split_reads_defined = 1
                split_reads = record.info["SR"]
            # create correct size
            if record.info["SVTYPE"] == "INS":
                if inserted_sequence == "N":
                    size = 5000
                else:
                    size = len(inserted_sequence)
            else:
                size = record.stop - record.pos
            # format output line
            output_line_list = [record.chrom, str(record.pos),
                                str(record.stop), ".", "-1", "-1",
                                str(record.id), "-1", ".", ".",
                                record.info["SVTYPE"], inserted_sequence, tool, str(read_pairs_defined),
                                str(read_pairs), str(split_reads_defined), str(split_reads),
                                str(size)]
            output_line = "\t".join(output_line_list) + "\n"
            # write output line
            bedpe.write(output_line)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # convert vcf file
    process_vcf(args.input_vcf, args.output_bedpe, args.tool)

if __name__ == "__main__":
    main()