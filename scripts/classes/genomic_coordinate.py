"""
Class for depicting genomic coordinate
"""

class GenomicCoordinate(object):

    def __init__(self, chromosome: str, position: int):
        """
        Genomic coordinate are 1-based
        """
        self.chromosome = chromosome
        self.position = position