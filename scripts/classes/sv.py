"""
Class to represent structural variant
"""

class SV(object):

    def __init__(self, chromosomeA: str, startA: int, endA: int, sv_type: str,
                 read_pairs: str = '0', split_reads: str = '0', chromosomeB: str = ".",
                 startB: int = -1, endB: int = -1, identifier: str = ".",
                 qual: int = -1, strandA: str = ".", strandB: str = ".",
                 tool: str = '', read_pairs_defined: int = 0,
                 split_reads_defined: int = 0, delly_low_qual_defined: int = 0,
                 delly_low_qual: int = 0, manta_ref_read_pairs_defined: int = 0,
                 manta_ref_read_pairs: int = 0,
                 manta_ref_split_reads_defined: int = 0,
                 manta_ref_split_reads: int = 0, inserted_sequence: str = "."):
        # coordinates are 1-end and inclusive
        self.chromosomeA = chromosomeA
        self.startA = startA
        self.endA = endA
        self.sv_type = sv_type
        self.chromosomeB = chromosomeB
        self.startB = startB
        self.endB = endB
        self.identifier = identifier
        self.qual = qual
        self.strandA = strandA
        self.strandB = strandB
        self.tool = tool
        self.read_pairs_defined = int(read_pairs_defined)
        self.read_pairs = read_pairs
        self.split_reads_defined = int(split_reads_defined)
        self.split_reads = split_reads
        self.delly_low_qual_defined = int(delly_low_qual_defined)
        self.delly_low_qual = int(delly_low_qual)
        self.manta_ref_read_pairs_defined = int(manta_ref_read_pairs_defined)
        self.manta_ref_read_pairs = int(manta_ref_read_pairs)
        self.manta_ref_split_reads_defined = int(manta_ref_split_reads_defined)
        self.manta_ref_split_reads = int(manta_ref_split_reads)
        self.inserted_sequence = inserted_sequence