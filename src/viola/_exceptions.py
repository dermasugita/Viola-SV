class TableNotFoundError(KeyError):
    pass

class InfoNotFoundError(KeyError):
    pass

class ContigNotFoundError(KeyError):
    pass

class VcfParseError(Exception):
    pass