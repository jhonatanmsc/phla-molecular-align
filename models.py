# import pdb

class Molecule():
    '''### Molecule Object
        #### attributes:
        - dbname: data base name
        - name: molecule name
        - seq: sequence of aminoacids
    '''

    def __init__(self, dbname: str, name: str, seq: str):
        self.dbname = dbname
        self.name = name
        self.seq = seq

    def __str__(self):
        return f'<Molecule {self.dbname} {self.name}>'
