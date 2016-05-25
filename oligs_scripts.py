# import os
import sys

from aux import *


def _examine_oligs_set(args):
    """
    Validate assemblies oligssetswise and show redundant oligs - whether they are main or mutant (if any).
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
    if len(used_oligs) < len(oligs) + len(mutant_oligs):
        red_oligs = filter(lambda x: x not in used_oligs, oligs)
        red_mutant_oligs = filter(lambda x: x not in used_oligs, mutant_oligs)
        if len(red_oligs) > 0:
            print 'Redundant main oligs:\n{}'.format('\n'.join(red_oligs))
        if len(red_mutant_oligs) > 0:
            print 'Redundant mutant oligs:\n{}'.format('\n'.join(red_mutant_oligs))
    used_old_oligs = filter(lambda x: x in used_oligs, old_oligs)
    used_old_mutant_oligs = filter(lambda x: x in used_oligs, old_mutant_oligs)
    print 'Number of used old oligs: {}'.format(len(used_old_oligs))
    print 'Number of used old mutant oligs: {}'.format(len(used_old_mutant_oligs))

    print '\n'


def main(args):
    option = args[0]
    if option == '--examine_set':
        return _examine_oligs_set(args[1:])
    raise KeyError('Unknown option')


if __name__ == '__main__':
    main(sys.argv[1:])
