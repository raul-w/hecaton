"""
Class to represent breakends of VCF files
"""

from classes.genomic_coordinate import GenomicCoordinate
from pysam import VariantRecord


class Breakend(object):

    def __init__(self, identifier: int, coordinate: GenomicCoordinate,
                 mate_coordinate: GenomicCoordinate,
                 mate_direction: str, mate_position: str,
                 variant_record: VariantRecord, mate_id: int = None,
                 inserted_sequence: str = "."):
        directions = ["L", "R"]
        if mate_direction not in directions:
            raise ValueError("Unknown orientation")
        if mate_position not in directions:
            raise ValueError("Unknown positions")
        self.identifier = identifier
        self.coordinate = coordinate
        self.inserted_sequence = inserted_sequence
        self.mate_id = mate_id
        self.mate_coordinate = mate_coordinate
        self.mate_direction = mate_direction
        self.mate_position = mate_position
        self.variant_record = variant_record
