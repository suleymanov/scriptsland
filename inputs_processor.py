import os
import sys

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

read_fasta = lambda x: list(SeqIO.parse(x, 'fasta'))

def _write_fasta(fn, records):
	with open(fn, 'w') as f:
		SeqIO.write(records, f, 'fasta')


def _check(args):
	excel_fn = args[0]
	fasta_fn = args[1]
	print 'Checking CD47 chain consistence in files {} and {}\n'.format(excel_fn, fasta_fn)

	df = pd.read_excel(excel_fn, sheetname=0)
	records = read_fasta(fasta_fn)

	names = list(df.iloc[:, 0])
	vl = list(df.iloc[:, 1])
	vh = list(df.iloc[:, 2])

	for i, name in enumerate(names):
		print 'Variant name: {}'.format(name)
		h_rec_name = name + ':H'
		h_rec = filter(lambda x: x.description == h_rec_name, records)
		assert len(h_rec) == 1
		print 'Heavy chains equal: {}'.format(vh[i] == str(h_rec[0].seq))
		l_rec_name = name + ':L'
		l_rec = filter(lambda x: x.description == l_rec_name, records)
		assert len(l_rec) == 1
		print 'Light chains equal: {}\n'.format(vl[i] == str(l_rec[0].seq))


def _check2(args):
	excel_fn = args[0]
	by_chain_pn = args[1]
	data = map(
		lambda x: {'name': x.split(os.sep)[-1].split('.')[0], 
		'sequence': str(read_fasta(by_chain_pn + os.sep + x)[0].seq)}, filter(
			lambda x: x.split('.')[-1] == 'fasta', 
			os.listdir(by_chain_pn)
		)
	)
	df = pd.read_excel(excel_fn, sheetname=0)
	names = list(df.iloc[:, 0])
	vl = list(df.iloc[:, 1])
	vh = list(df.iloc[:, 2])

	for i, name in enumerate(names):
		items = filter(lambda x: name in x['name'], data)
		assert len(items) == 2
		assert vl[i] in map(lambda x: x['sequence'], items)
		assert vh[i] in map(lambda x: x['sequence'], items)


def _create_new(args):
	excel_fn = args[0]
	by_chain_pn = args[1]
	if os.path.exists(by_chain_pn):
		os.rmdir(by_chain_pn)
	os.mkdir(by_chain_pn)
	df = pd.read_excel(excel_fn, sheetname=0)
	names = list(df.iloc[:, 0])
	vl = list(df.iloc[:, 1])
	vh = list(df.iloc[:, 2])
	for i, name in enumerate(names):
		h_name = name + ':H'
		h_fn_name = by_chain_pn + os.sep + name + '_vh.fasta'
		rec = SeqRecord(Seq(vh[i]), name=h_name, id=h_name, description='')
		_write_fasta(h_fn_name, [rec])

		l_name = name + ':L'
		l_fn_name = by_chain_pn + os.sep + name + '_vl.fasta'
		rec = SeqRecord(Seq(vl[i]), name=l_name, id=l_name, description='')
		_write_fasta(l_fn_name, [rec])


def main(args):
	option = args[0]
	if option == '--check':
		return _check(args[1:])
	if option == '--check2':
		return _check2(args[1:])
	if option == '--createnew':
		return _create_new(args[1:])
	raise ValueError('Wrong option.')


if __name__ == '__main__':
	main(sys.argv[1:])
