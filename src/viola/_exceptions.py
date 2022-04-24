class TableNotFoundError(KeyError):
    pass

class InfoNotFoundError(KeyError):
    pass
class SVIDNotFoundError(KeyError):
    pass

class ContigNotFoundError(KeyError):
    pass
class IllegalArgumentError(ValueError):
    pass
class DuplicatedPatientIDError(Exception):
    pass