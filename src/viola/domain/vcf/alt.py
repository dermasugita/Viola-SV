from pydantic import BaseModel

class Alt(BaseModel):
    _id: str
    description: str