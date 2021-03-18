class Indexer(object):
    @property
    def idx(self):
        return SvIdIndexer("idx", self)

class RootIndexer:
    def __init__(self, name, obj):
        self.name = name
        self.obj = obj

class SvIdIndexer(RootIndexer):
    def __getitem__(self, value):
        if type(value) is str:
            value = [value]
        return self.obj.filter_by_id(value)
    