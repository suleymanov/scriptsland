import sys

from utils import *


def _enlarge_last_oligs(args):
    """
    Correct set's latter oligs to improve amplification quality.
    :param args: list of str
    :return: None
    """
    min_len, max_len = 55, 59
    by_oligs_pn = args[0]
    oligs_fn = args[1]
    mutant_oligs_fn = args[2]
    variants_fn = args[3]
    set_pn = os.path.abspath(by_oligs_pn)

    print 'Processing oligs sets path: {}'.format(set_pn)
    oligs_sets = map(lambda x: SourceSet(x, read_fasta(x)), find_all_fasta(by_oligs_pn))
    oligs = read_fasta(oligs_fn)
    mut_oligs = read_fasta(mutant_oligs_fn)
    variants = map(get_rec_seq, read_fasta(variants_fn))

    for olig_set in oligs_sets:
        components = olig_set.components
        _, na_chain = assemble_chain(components)
        if len(components[-2]) < min_len:
            comp_str, comp_name = get_rec_seq(components[-2]), components[-2].id
            ind = na_chain.find(comp_str)
            if ind < 0:
                raise ValueError('Cannot find olig in sequence.')
            new_comp_str = na_chain[ind:ind + max_len]
            components[-2] = SeqRecord(Seq(new_comp_str), id=comp_name, name=comp_name, description='')
            if comp_str in map(get_rec_seq, oligs) and new_comp_str not in map(get_rec_seq, oligs):
                oligs.append(components[-2])
                print '\tAppended to main oligs.'
            if comp_str in map(get_rec_seq, mut_oligs) and new_comp_str not in map(get_rec_seq, mut_oligs):
                mut_oligs.append(components[-2])
                print '\tAppended to mutant oligs.'
            if assemble_chain(components)[0] not in variants:
                raise ValueError('Could not assemble chain with new oligs.')
        if len(components[-1]) < min_len:
            comp_str, comp_name = get_rec_seq(components[-1]), components[-1].id
            ind = get_rev_compl2(na_chain).find(comp_str)
            if ind < 0:
                raise ValueError('Cannot find olig in sequence.')
            new_comp_str = get_rev_compl2(na_chain)[ind:ind + max_len]
            components[-1] = SeqRecord(Seq(new_comp_str), id=comp_name, name=comp_name, description='')
            if comp_str in map(get_rec_seq, oligs) and new_comp_str not in map(get_rec_seq, oligs):
                oligs.append(components[-1])
                print '\tAppended to main oligs.'
            if comp_str in map(get_rec_seq, mut_oligs) and new_comp_str not in map(get_rec_seq, mut_oligs):
                mut_oligs.append(components[-1])
                print '\tAppended to mutant oligs.'
            if assemble_chain(components)[0] not in variants:
                raise ValueError('Could not assemble chain with new oligs.')
        write_fasta(olig_set.source_name, components)
    write_fasta(oligs_fn, oligs)
    write_fasta(mutant_oligs_fn, mut_oligs)
    print '\n'


def main(args):
    option = args[0]
    if option == '--enlarge':
        return _enlarge_last_oligs(args[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
