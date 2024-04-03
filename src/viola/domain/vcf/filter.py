from pydantic import BaseModel

class Filter(BaseModel):
    _id: str
    description: str