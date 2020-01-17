'''
    Script python for alignment of molecules
    ===================================================

    Author: Jhonatanmsc
    Github: https://github.com/jhonatanmsc/phla-molecular-align
'''


def load_molecules(filename: str):
    '''load molecules from fasta files
    :param filename: Your file name
    :returns: dict with molecules
    :raises: None
    '''
    molecules: list
    with open(filename) as file:
        lines = file.read()
        molecules_str = lines.slipt('>')
        for molecule in molecules_str:
            molecule = molecule.split('\n', 1)
            molecules.append({molecule[0]: molecule[1]})

    return molecules


def align_molecules(mol1: str, mol2: str):
    '''Do alignment of two molecules
    :params mol1, mol2: A name of your molecule
    :returns: dict with alignment result and indentity
    :raises: None
    '''
    
    return {}


def save_in_file(filename: str, alignment: dict):
    '''save result of alignment in file.txt
    :param filename: Your file name
    :param alignment: A dict with result of alignments of molecules
    :returns: None
    :raises: None
    '''

    with open(filename, 'w') as file:
        for molname in alignment:
            file.write(alignment[molname]['result'])

def main():
    molecules = {
        'p3d': load_molecules('molecules/phla3d.fasta'),
        'imgt': load_molecules('molecules/imgthla.fasta')
    }

    aligned = {}
    aligned_err = {}
    alignment = {}

    for molname in molecules['p3d']:
        mol1 = molecules['p3d'][molname]
        mol2 = molecules['imgt'][molname]
        alignment = align_molecules(mol1, mol2)

        if alignment['identity'] < 1:
            aligned_err[molname] = alignment
        else:
            aligned = alignment

    save_in_file('alignment.txt', aligned)
    save_in_file('alignment_err.txt', aligned_err)


if __name__ == '__main__':
    main()