from __future__ import annotations
from pydantic import BaseModel
from enum import Enum
from viola.domain.vcf.data_types import DataType
from viola.domain.vcf.number import Number



class Info(BaseModel):
    _id: str
    number: Number
    type: DataType
    description: str
    source: str | None
    version: str | None
    # Other fields in INFO column