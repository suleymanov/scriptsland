# import os
import sys

from itertools import izip

from utils import *


def _examine_oligs_set(args):
    """
    Validate assemblies oligssetswise, show redundant oligs and exclude them.
    :param args: list of str
    :return: None
    """
    by_oligs_pn = args[0]
    oligs_fn = args[1]
    mutant_oligs_fn = args[2]
    variants_fn = args[3]
    old_oligs_fn = args[4]
    old_mutant_oligs_fn = args[5]
    set_pn = os.path.abspath(by_oligs_pn)

    oligs_sets = map(lambda x: SourceSet(source_name=x, components=read_fasta(x)), find_all_fasta(by_oligs_pn))
    oligs = map(get_rec_seq, read_fasta(oligs_fn))
    mutant_oligs = map(get_rec_seq, read_fasta(mutant_oligs_fn))
    variants = map(get_rec_seq, read_fasta(variants_fn))
    old_oligs = map(get_rec_seq, read_fasta(old_oligs_fn))
    old_mutant_oligs = map(get_rec_seq, read_fasta(old_mutant_oligs_fn))

    print 'Set directory: {}'.format(set_pn)
    print 'Number of main oligs: {}'.format(len(oligs))
    print 'Number of unique main oligs: {}'.format(len(set(oligs)))
    print 'Number of mutant oligs: {}'.format(len(mutant_oligs))
    print 'Number of unique mutant oligs: {}'.format(len(set(mutant_oligs)))
    print 'Number of variants: {}'.format(len(variants))
    intersection = set(oligs).intersection(mutant_oligs)
    print 'Intersection length: {}'.format(len(intersection))
    if len(intersection) > 0:
        print '\n'.join(intersection)
    print 'Number of old oligs used: {}'.format(len(set(old_oligs).intersection(oligs)))
    print 'Number of old mutant oligs used: {}'.format(len(set(old_mutant_oligs).intersection(mutant_oligs)))

    assemblies = map(assemble_chain, map(get_attr('components'), oligs_sets))
    aa_assemblies, na_assemblies = map(lambda x: x[0], assemblies), map(lambda x: x[1], assemblies)
    assert len(aa_assemblies) == len(variants) and all(s in aa_assemblies for s in variants) and all(s in variants for s in aa_assemblies)

    used_oligs = list(set(reduce(lambda res, item: res + map(get_rec_seq, item.components), oligs_sets, [])))
    print 'Number of used oligs: {}'.format(len(used_oligs))
    mutant_oligs = filter(lambda x: x not in oligs, mutant_oligs)
    if len(used_oligs) < len(oligs) + len(mutant_oligs):
        red_oligs = filter(lambda x: x not in used_oligs, oligs)
        red_mutant_oligs = filter(lambda x: x not in used_oligs, mutant_oligs)
        if len(red_oligs) > 0:
            print 'Redundant main oligs:\n{}'.format('\n'.join(red_oligs))
            write_fasta(oligs_fn, filter(lambda x: get_rec_seq(x) not in red_oligs, read_fasta(oligs_fn)))
        if len(red_mutant_oligs) > 0:
            print 'Redundant mutant oligs:\n{}'.format('\n'.join(red_mutant_oligs))
            write_fasta(mutant_oligs_fn, filter(lambda x: get_rec_seq(x) not in red_mutant_oligs, read_fasta(mutant_oligs_fn)))
    used_old_oligs = filter(lambda x: x in used_oligs, old_oligs)
    used_old_mutant_oligs = filter(lambda x: x in used_oligs, old_mutant_oligs)
    print 'Number of used old oligs: {}'.format(len(used_old_oligs))
    print 'Number of used old mutant oligs: {}'.format(len(used_old_mutant_oligs))

    print '\n'


def _final_set_improvement(args):
    """
    - minimize number of main/mutant oligs
    - validate if all assemblies are correct
    :param args: list of str
    :return: None
    """
    oligs_fn = args[0]
    mut_fn = args[1]
    by_oligs_pn = args[2]
    variants_fn = args[3]
    set_pn = os.path.abspath(by_oligs_pn)

    oligs = map(get_rec_seq, read_fasta(oligs_fn))
    mut_oligs = map(get_rec_seq, read_fasta(mut_fn))
    oligs_sets = map(lambda x: SourceSet(x, read_fasta(x)), find_all_fasta(by_oligs_pn))
    variants = map(get_rec_seq, read_fasta(variants_fn))

    print 'Set directory: {}'.format(set_pn)
    print 'Original number of main oligs: {}'.format(len(oligs))
    print 'Original number of mutant oligs: {}'.format(len(mut_oligs))

    mut_oligs = filter(lambda x: x not in oligs, mut_oligs)
    used_oligs = list(set(reduce(lambda res, item: res + map(get_rec_seq, item.components), oligs_sets, [])))
    assemblies = map(lambda x: x[0], map(assemble_chain, map(get_attr('components'), oligs_sets)))
    assert len(oligs) == len(set(oligs))
    assert len(mut_oligs) == len(set(mut_oligs))
    assert len(oligs) + len(mut_oligs) == len(used_oligs)
    assert len(assemblies) == len(variants)
    assert all(s in variants for s in assemblies)
    assert all(s in assemblies for s in variants)
    assert len(assemblies) == len(set(assemblies))


def main(args):
    option = args[0]
    if option == '--examine_set':
        return _examine_oligs_set(args[1:])
    if option == '--improve':
        return _final_set_improvement(args[1:])
    raise KeyError('Unknown option')


if __name__ == '__main__':
    main(sys.argv[1:])
