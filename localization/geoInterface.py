""" Anchor and target class"""


class Anchor:
    def __init__(self, ID, loc):
        self.loc = loc
        self.ID = str(ID)

    def __str__(self):
        return 'Anchor ' + self.ID + ' @ ' + self.loc.__str__()


class Target:
    def __init__(self, ID):
        self.loc = None
        self.ID = str(ID)
        self.measures = []

    def __str__(self):
        if self.loc is None:
            return 'Target ' + self.ID
        else:
            return 'Target ' + self.ID + ' @ Real Location:' + self.loc.__str__()

    def add_measure(self, a, d):
        self.measures.append((a, d))
