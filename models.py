# import pdb


class Resi():
    '''### Resi Object
        #### attributes:
        - pos: aminoacid position
        - name: aminoacid name abreviated
    '''

    def __init__(self, pos: int, name: str):
        self.pos = pos
        self.name = name

    def __str__(self):
        return f'<Resi {self.pos} {self.name}'

    def __eq__(self, other):
        return self.pos == other.pos and self.name == other.name


class Molecule():
    '''### Molecule Object
        #### attributes:
        - dbname: data base name
        - name: molecule name
        - seq: sequence of aminoacids
    '''

    def __init__(self, dbname: str, name: str, seq: list):
        self.dbname = dbname
        self.name = name
        self.seq = seq

    def __str__(self):
        return f'<Molecule {self.dbname} {self.name}>'

    def resi_name(self, pos: int):
        "#### receive position and returns molecule"
        for resi in self.seq:
            if resi.pos == pos:
                return resi
        return None

    def all_pos(self):
        "#### returns all molecule positions"
        list_pos = []
        for mol in self.seq:
            list_pos.append(mol.pos)
        return list_pos
