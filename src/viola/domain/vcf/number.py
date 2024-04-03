from pydantic import BaseModel

class Number(BaseModel):
    value: str
    ALLOWED_CHARACTERS = "ARG."
    def __init__(self, value: str):
        if not isinstance(value, (int, str)) or (isinstance(value, str) and value not in self.ALLOWED_CHARACTERS):
            message = f"Value must be an integer or one of 'A', 'R', 'G', or '.'"
            raise ValueError(message)
        self.value = value

    def __repr__(self):
        return f"Number({self.value})"
    
    def __str__(self):
        return self.value