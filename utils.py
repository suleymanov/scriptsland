import os

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import unambiguous_dna

from collections import namedtuple
from operator import attrgetter, itemgetter
from simplejson import dumps, loads


Component = namedtuple('Component', ['plate_name', 'pos', 'name', 'seq', 'vol'])
SourceSet = namedtuple('SourceSet', ['source_name', 'components'])

get_attr = lambda key: attrgetter(key)
get_item = lambda key: itemgetter(key)
get_short_fn = lambda x: x.split(os.sep)[-1].split('.')[0]
create_seq = lambda x: ''.join(np.random.choice(list(unambiguous_dna.letters), x, replace=True))
create_seq2 = lambda x: Seq(x, generic_dna)
create_record = lambda i, s: SeqRecord(seq=Seq(s), name='Seq_{}'.format(i), id='Seq_{}'.format(i), description='')
create_olig_record = lambda i, s, num_well: SeqRecord(seq=Seq(s), name='Olig_{}'.format(i), id='Olig_{}'.format(i), description='{}'.format(num_well))
create_olig_record2 = lambda s, beg, end: SeqRecord(seq=Seq(s), name='olig_[{}, {})'.format(beg, end), id='olig_[{}, {})'.format(beg, end), description='')
read_fasta = lambda x: list(SeqIO.parse(x, 'fasta'))
find_all_fasta = lambda x: [x + os.sep + fn for fn in os.listdir(x) if fn.split('.')[-1] == 'fasta']
get_rec_seq = lambda x: str(x.seq)  # Bio.SeqRecord.SeqRecord -> str
get_rev_compl = lambda x: str(x.seq.reverse_complement())  # Bio.SeqRecord.SeqRecord -> str
get_rev_compl2 = lambda x: str(Seq(x).reverse_complement())  # str -> str
remove_gaps = lambda x: SeqRecord(seq=Seq(str(x.seq).replace('-', '')), name=x.description, id=x.description, description='')  # Bio.SeqRecord.SeqRecord -> BioSeqRecord.SeqRecord


def assemble_chain(olig_records):
    get_common_inds = lambda s1, s2: [i for i in xrange(min(len(s1), len(s2)))
                        if s2[:(i + 1)] == s1[(-i - 1):]]
    get_next_index = lambda s1, s2: (get_common_inds(s1, s2)[-1] + 1
                        if len(get_common_inds(s1, s2)) > 0 else 0)
    na_chain = ''
    for i in xrange(len(olig_records)):
        next_olig = (get_rec_seq(olig_records[i])
                        if i % 2 == 0 else get_rev_compl(olig_records[i]))
        na_chain += next_olig[get_next_index(na_chain, next_olig):]
    return str(Seq(na_chain).translate()), na_chain


def read_json(fn):
    with open(fn) as f:
        return loads(f.read())


def write_json(fn, data):
    with open(fn, 'w') as f:
        f.write(dumps(data, indent=4))


def write_fasta(fn, data):
    with open(fn, 'w') as f:
        SeqIO.write(data, f, 'fasta')


def write_file(fn, data):
    with open(fn, 'w') as f:
        f.write(data)


def read_file(fn):
    with open(fn) as f:
        return f.read()


def lcsubstr(str1, str2):
    """
    Returns longest common substring of two strings.
    """
    shape = (len(str1), len(str2))
    aux = np.zeros(shape)
    z = 0
    res = []
    get_substr = lambda s, i, z: s[(i-z+1):(i+1)]
    for i in xrange(shape[0]):
        for j in xrange(shape[1]):
            if str1[i] == str2[j]:
                if i == 0 or j == 0:
                    aux[i, j] = 1
                else:
                    aux[i, j] = aux[i-1, j-1] + 1
                if aux[i, j] > z:
                    z = int(aux[i, j])
                    res = [get_substr(str1, i, z)]
                elif aux[i, j] == z:
                    res.append(get_substr(str1, i, z))
            else:
                aux[i, j] = 0
    return res[0] if len(res) > 0 else ''