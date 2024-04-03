from pydantic import BaseModel
from viola.domain.vcf.data_types import DataType
from viola.domain.vcf.number import Number

class Format(BaseModel):
    _id: str
    number: Number
    type: DataType
    description: str