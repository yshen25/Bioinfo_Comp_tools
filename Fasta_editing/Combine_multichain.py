from Bio import SeqIO

def rename(full_name:str) -> str:
    return full_name.replace('*', '').replace(':','_')

def filename(A_name, B_name):
    A_new = rename(A_name)
    B_new = rename(B_name)

    return f"{A_new}_{B_new}.faa"

def AB_combine(A_fasta, A_allele_file, B_fasta, B_allele_file):
    """
    Read in two fasta files containing A chain sequences and B chain sequences,
    combining both chains into single fasta files according to A and B allele files
    """
    A_record = SeqIO.parse(A_fasta, 'fasta')
    B_record = SeqIO.parse(B_fasta, 'fasta')

    with open(A_allele_file) as fh:
        A_list = [line.strip() for line in fh]

    with open(B_allele_file) as fh:
        B_list = [line.strip() for line in fh]

    for A_allele in A_list:
        for B_allele in B_list:
            A_rec = A_record[A_allele]
            B_rec = B_record[B_allele]

            SeqIO.write([A_rec, B_rec], filename(A_allele, B_allele), 'fasta')

    return

if __name__ == '__main__':
    AB_combine()