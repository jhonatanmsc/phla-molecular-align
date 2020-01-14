'''
    Script python for alignment of molecules
    ===================================================

    Author: Jhonatanmsc
    Github: https://github.com/jhonatanmsc/phla-molecular-align
'''

# load molecules from fasta files
# returns dict with molecules
def load_molecules(filename):
    return {}

# Do alignment of two molecules
# returns dict with alignment result and identity
def align_molecules(mol1, mol2):
    return {}

# save result of alignment in file.txt
def save_in_file(filename, alignment):
    file = open(filename, 'w')
    for molname in alignment:
        file.write(alignment[molname]['result'])
    file.close()
    return file

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