'''
    Script python for alignment of molecules
    ===================================================

    Author: Jhonatanmsc
    Github: https://github.com/jhonatanmsc/phla-molecular-align
'''
from models import Molecule
from Bio import pairwise2
from Bio.pairwise2 import align
from Bio.SubsMat.MatrixInfo import blosum62
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

            ignorable_alleles = ['N', 'L', 'Q', 'S', 'A', 'C']

            name = mol_str[0].split(' ')
            if len(name) > 2:
                name = name[1]
            elif len(name) > 1:
                name = name[0]
            else:
                name = name[0]

            if any(name.endswith(allele) for allele in ignorable_alleles):
                print('ignorated ', name)
                continue
            
            if name:
                if name.find(':') > 2:
                    name = ':'.join(name.split(':', 2)[:2])

                if any(subname in molecules for subname in name):
                    continue

            else:
                name = 'None'

            mol = Molecule(dbname=dbname, name=name, seq=mol_str[1])
            molecules[name] = mol

    return molecules


def format_alignment(mol1: Molecule, mol2: Molecule):
    '''### Do alignment of two molecules
        #### params:
        - mol1, mol2: You includer molecule to align

        *returns* -> Molecule with alignment result and indentity
    '''
    alignment = align.localds(mol1.seq, mol2.seq,  blosum62, -12, -4)
    
    alignment_formated = pairwise2.format_alignment(*alignment[0])
    alignment_formated = alignment_formated.split('\n')
    
    header = ''
    seq_mol1 = alignment_formated[0]
    seq_mol2 = alignment_formated[2]
    result_raw = alignment_formated[1]
    identity = alignment[0][-1]

    body_mol1 = ''
    body_mol2 = ''
    result = ''
    body = ''
    count = 0
    errors = 0

    for i in range(len(seq_mol1)):
        body_mol1 += seq_mol1[i]
        body_mol2 += seq_mol2[i]
        result += result_raw[i]

        if not seq_mol1[i].isnumeric() and seq_mol1[i] != ' ':
            if seq_mol1[i] == seq_mol2[i]:
                count += 1
            else:
                errors += 1

        if (i+1) % 60 == 0:
            body += f"{body_mol1}\n{result}\n{body_mol2}\n\n"
            result = ''
            body_mol1 = ''
            body_mol2 = ''

    identity = count/(count + errors)
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
        if molname in molecules['p3d'] and molname in molecules['imgt']:
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
