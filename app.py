'''
    Script python for alignment of molecules
    ===================================================

    Author: Jhonatanmsc
    Github: https://github.com/jhonatanmsc/phla-molecular-align
'''
from models import Molecule, Resi
import re
# import pdb


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
            pre_gsh = mol_str[1][:mol_str[1].find('GSH')]
            pos_gsh = mol_str[1][mol_str[1].find('GSH'):]

            seq = []
            for i in range(len(pre_gsh)):
                resi = Resi(pos=i-len(pre_gsh), name=pre_gsh[i])
                seq.append(resi)

            for i in range(len(pos_gsh)):
                resi = Resi(pos=i, name=pos_gsh[i])
                seq.append(resi)
            name = re.search(r'[A-Z]\*[\d]{2}\:[\d]{2}', mol_str[0])
            if name:
                name = name.group()
            else:
                name = ''
            mol = Molecule(dbname=dbname, name=name, seq=seq)
            molecules[name] = mol

    return molecules


def align_molecules(mol1: Molecule, mol2: Molecule):
    '''### Do alignment of two molecules
        #### params:
        - mol1, mol2: Your molecule to align
        
        *returns* -> Molecule with alignment result and indentity
    '''
    header = ''
    seq_mol1 = ''
    seq_mol2 = ''
    result = ''
    body = ''
    lines = 1
    count = 0
    positions = mol1.all_pos() + list(set(mol2.all_pos()) - set(mol1.all_pos()))
    positions.sort()

    for pos in positions:
        mol1_resi = mol1.resi_name(pos)
        mol2_resi = mol2.resi_name(pos)
        if mol1_resi and mol2_resi:
            if mol1_resi == mol2_resi:
                result += '*'
                count += 1
            else:
                result += '_'
            seq_mol1 += mol1_resi.name
            seq_mol2 += mol2_resi.name
        elif mol1_resi:
            seq_mol1 += mol1_resi.name
            seq_mol2 += ' '
            result += '_'
        elif mol2_resi:
            seq_mol1 += ' '
            seq_mol2 += mol2_resi.name
            result += '_'
        else:
            result += ' '

        if lines % 60 == 0:
            nw_lines = len(str(lines-59))
            nq_lines = 8-nw_lines
            n_lines = (' ' * nq_lines) + str(lines-59)
            body += f"{n_lines} {seq_mol1}\n"
            body += f"{n_lines} {seq_mol2}\n"
            body += f"{'        '} {result}\n"

            result = ''
            seq_mol1 = ''
            seq_mol2 = ''
        lines += 1

    header = "< %s - %s | %s | %d%%\n" % (mol1.dbname,
                                          mol2.dbname, mol1.name, (count / len(mol1.seq)) * 100)
    text = header + body
    
    return {
        'text': text,
        'identity': (count / len(mol1.seq))
    }


def main():
    molecules = {
        'p3d': load_molecules('molecules/phla3d.fasta', 'p3d'),
        'imgt': load_molecules('molecules/imgthla.fasta', 'imgt')
    }

    alignments = []
    for molname in molecules['p3d']:
        mol1 = molecules['p3d'][molname]
        mol2 = molecules['imgt'][molname]
        result = align_molecules(mol1, mol2)
        alignments.append(result)

    alignment_ok = ''
    alignment_err = ''
    for alignment in alignments:
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
