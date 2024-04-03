from pydantic import BaseModel

class Contig(BaseModel):
    _id: str
    length: int | None
    assembly: str | None
    md5: str | None
    species: str | None
    taxonomy: str | None
