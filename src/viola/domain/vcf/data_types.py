from __future__ import annotations
from enum import Enum

class DataType(str, Enum):
    """
    reference: https://samtools.github.io/hts-specs/VCFv4.3.pdf
    """
    Integer = "Integer"
    Float = "Float"
    Flag = "Flag"
    String = "String"
    Character = "Character"
    def __repr__(self):
        return f"DataType.{self.value}"
    def __str__(self):
        return self.value
    
    def __eq__(self, other: DataType):
        return self.value == other.value
    
    def __ne__(self, other: DataType):
        return self.value != other.value