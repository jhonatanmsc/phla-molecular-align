'''
    Script python for alignment of molecules
    ===================================================

    Author: Jhonatanmsc
    Github: https://github.com/jhonatanmsc/phla-molecular-align
'''
from models import Molecule
from Bio import pairwise2
from Bio.pairwise2 import align
import pdb


def load_molecules(filename: str, dbname='#'):
    '''### load molecules from fasta file
        #### params:
        - filename: Your file name
        - dbname: Name of your data base

        *returns* -> dict with molecules
    '''
    molecules = {}
    with open(filename) as file:
        lines = file.read()
        molecules_str = lines.split('>')[1::]

        for molecule_str in molecules_str:
            mol_str = molecule_str.split('\n', 1)
            mol_str[1] = mol_str[1].replace('\n', '')

            # ignore null alleles
            if mol_str[0].find('N') == 0:
                continue

            # ignore lower alleles
            if mol_str[0].find('L') == 0:
                continue

            # ignore questionable alleles
            if mol_str[0].find('Q') == 0:
                continue

            # ignore secreted
            if mol_str[0].find('S') == 0:
                continue

            # ignore aberrant
            if mol_str[0].find('A') == 0 and mol_str[0].find('C') == 0:
                continue

            name = mol_str[0].split(
                ' ')[1] if ' ' in mol_str[0] else mol_str[0]

            if name:
                if name.find(':') > 2:
                    name = name.split(':')[0] + ':' + name.split(':')[1]

                if name in molecules:
                    continue

            else:
                name = ''

            mol = Molecule(dbname=dbname, name=name, seq=mol_str[1])
            molecules[name] = mol

    return molecules


def format_alignment(mol1: Molecule, mol2: Molecule):
    '''### Do alignment of two molecules
        #### params:
        - mol1, mol2: Youincluder molecule to align

        *returns* -> Molecule with alignment result and indentity
    '''
    alignment = align.globalxx(mol1.seq, mol2.seq)[0]
    header = ''
    seq_mol1 = ''
    seq_mol2 = ''
    result = ''
    body = ''
    count = 0

    for i in range(alignment[4]):
        seq_mol1 += alignment[0][i]
        seq_mol2 += alignment[1][i]

        if alignment[0][i] == alignment[1][i]:
            count += 1
            result += '|'
        else:
            result += ' '

        if (i+1) % 60 == 0:
            body += f"{seq_mol1}\n{result}\n{seq_mol2}\n\n"
            result = ''
            seq_mol1 = ''
            seq_mol2 = ''

    identity = len(mol1.seq)/count
    pdb.set_trace()
    header = "< %s - %s | %s | %.1f%%\n" % (mol1.dbname,
                                            mol2.dbname, mol1.name, identity * 100)
    text = header + body

    return {
        'text': text,
        'identity': identity
    }


def main():
    molecules = {
        'p3d': load_molecules('molecules/phla3d.fasta', 'p3d'),
        'imgt': load_molecules('molecules/imgthla.fasta', 'imgt'),
    }

    alignments_formated = []
    for molname in molecules['p3d']:
        alignments_formated.append(format_alignment(
            molecules['p3d'][molname], molecules['imgt'][molname]))

    alignment_ok = ''
    alignment_err = ''
    for alignment in alignments_formated:
        if alignment['identity'] == 1:
            alignment_ok += alignment['text']
        else:
            alignment_err += alignment['text']

    if alignment_err != '':
        with open('alignment-err.txt', 'w') as file:
            file.write(alignment_err)
    if alignment_ok != '':
        with open('alignment-ok.txt', 'w') as file:
            file.write(alignment_ok)


if __name__ == '__main__':
    main()
