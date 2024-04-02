from __future__ import annotations
from pydantic import BaseModel
from enum import Enum

class InfoNumber(BaseModel):
    value: str
    ALLOWED_CHARACTERS = "ARG."
    def __init__(self, value: str):
        if not isinstance(value, (int, str)) or (isinstance(value, str) and value not in self.ALLOWED_CHARACTERS):
            message = f"Value must be an integer or one of 'A', 'R', 'G', or '.'"
            raise ValueError(message)
        self.value = value

    def __repr__(self):
        return f"InfoNumber({self.value})"
    
    def __str__(self):
        return self.value

class InfoType(str, Enum):
    INTEGER = "Integer"
    FLOAT = "Float"
    FLAG = "Flag"
    STRING = "String"

    def __repr__(self):
        return f"InfoType.{self.value}"
    def __str__(self):
        return self.value


class Info(BaseModel):
    _id: str
    number: InfoNumber
    type: InfoType
    description: str
    source: str | None
    version: str | None
    # Other fields in INFO column