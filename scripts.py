import sys

from aux import *


def _enlarge_last_oligs(args):
    """
    Correct set's latter oligs to improve amplification quality.
    :param args: list of str
    :return: None
    """
    by_oligs_pn = args[0]
    oligs_fn = args[1]
    mutant_oligs_fn = args[2]
    variants_fn = args[3]

    oligs_sets = map(lambda x: SourceSet(x, read_fasta(x)), find_all_fasta(by_oligs_pn))
    oligs = map(get_rec_seq, read_fasta(oligs_fn))
    mut_oligs = map(get_rec_seq, read_fasta(mutant_oligs_fn))
    variants = map(get_rec_seq, read_fasta(variants_fn))

    for olig_set in oligs_sets:
        components = olig_set.components


def main(args):
    raise NotImplementedError()


if __name__ == '__main__':
    main(sys.argv[1:])
