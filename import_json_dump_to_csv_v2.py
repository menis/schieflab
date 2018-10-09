#!/Users/menis/anaconda/bin/python
from abutils.core import sequence
from abutils.core import pair
import pandas as pd
import numpy as np
import json
from pandas.io.json import json_normalize
import sys, os

pairs = []
with open(sys.argv[1]) as f:
    for line in f:
        d = json.loads(line.strip())
        # Bryan has a pair object
        seq = sequence.Sequence(d)
        seqs = [seq]
        abpair = pair.Pair(seqs)
        pairs.append(abpair)


# Built in json flatten from pandas does an OK job.  
#json_df = json_normalize(d)
#json_df.to_csv("dump.csv", sep=',')

# An example of a mole elaborate json flatten that we can repurpose
# From https://towardsdatascience.com/flattening-json-objects-in-python-f5343c794b10
def flatten_json(y):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


# Extracted from Bryan's code to fix output_schief_csv()

# FROM https://github.com/briney/abtools/blob/master/abtools/pipeline.py
def make_dir(d):
    '''
    Makes a directory, if it doesn't already exist.
    Args:
        d (str): Path to a directory.
    '''
    if not os.path.exists(d):
        os.makedirs(d)

# FROM https://github.com/briney/vaxtools/blob/master/vaxtools/utils/outputs.py
def schief_csv_output(pairs, output_file, sep=',', legacy_abstar=True):
    make_dir(os.path.dirname(output_file))
    header = _get_schief_output_header(sep)
    output = [header, ]
    for p in sorted(pairs, key=lambda x: _get_name(x)):
        name = _get_name(p)
        line = [name, ]
        #line += _get_pair_metadata(p)
        line += _schief_output_line(p.heavy, legacy_abstar)
        line += _schief_output_line(p.light, legacy_abstar)
        output.append(sep.join([str(l) for l in line]))
    open(output_file, 'w').write('\n'.join(output))


def _get_name(p):
    name = ''
    if p.heavy is not None:
        if 'seq_id' in p.heavy:
            name = p.heavy['seq_id']
    elif p.light is not None:
        if 'seq_id' in p.light:
            name = p.light['seq_id']
    return name if name != '' else p.name


def _get_pair_metadata(p):
    if p.heavy is not None:
        seq = p.heavy
    else:
        seq = p.light
    experiment = seq.get('experiment', '')
    group = seq.get('group', '')
    subject = seq.get('subject', '')
    timepoint = seq.get('timepoint', '')
    # experiment = seq['experiment'] if 'experiment' in seq else ''
    # group = seq['group'] if 'group' in seq else ''
    # subject = seq['subject'] if 'subject' in seq else ''
    # timepoint = seq['timepoint'] if 'timepoint' in seq else ''
    return [experiment, group, subject, timepoint]


def _get_fr_identity(seq, res='nt'):
    len_field = 'region_len_nt' if res == 'nt' else 'region_len_aa'
    mut_field = 'region_muts_nt' if res == 'nt' else 'region_muts_aa'
    regions = ['fr1', 'fr2', 'fr3']
    length = sum([seq[len_field][region] for region in regions])
    muts = sum([len(seq[mut_field][region]['muts']) for region in regions])
    return 100. * muts / length


def _get_schief_output_header(sep):
 '''
   fields = ['Sequence ID', 'Experiment', 'Group', 'Subject', 'Timepoint', 'VH gene', 'DH gene', 'JH gene', 'CDR3 length',
              'Junction AA', 'Junction NT seq', '% VH mutation (NT)', '% FR mutation (NT)',
              '% VH mutation (AA)', '% FR mutation (AA)', 'VH insertions', 'VH deletions',
              'VDJ AA seq', 'VDJ NT seq', 'Insertion count', 'Insertion lengths',
              'Deletion count', 'Deletion lengths', 'VDJ cysteine count', 'CDR3 cysteine count',
              'VL gene', 'JL gene', 'CDR3 length', 'Junction AA', 'Junction NT seq',
              '% VL mutation (NT)', '% FR mutation (NT)', '% VL mutation (AA)', '% FR mutation (AA)',
              'VL insertions', 'VL deletions', 'VDJ AA seq', 'VDJ NT seq',
              'Insertion count', 'Insertion lengths', 'Deletion count', 'Deletion lengths',
              'VDJ cysteine count', 'CDR3 cysteine count']
 '''
 fields = ['Sequence ID', 'VH gene', 'DH gene', 'JH gene',
              'CDR3 length',
              'Junction AA', 'Junction NT seq', '% VH mutation (NT)', '% FR mutation (NT)',
              '% VH mutation (AA)', '% FR mutation (AA)', 'VH insertions', 'VH deletions',
              'VDJ AA seq', 'VDJ NT seq', 'Insertion count', 'Insertion lengths',
              'Deletion count', 'Deletion lengths', 'VDJ cysteine count', 'CDR3 cysteine count',
              'VL gene', 'JL gene', 'CDR3 length', 'Junction AA', 'Junction NT seq',
              '% VL mutation (NT)', '% FR mutation (NT)', '% VL mutation (AA)', '% FR mutation (AA)',
              'VL insertions', 'VL deletions', 'VDJ AA seq', 'VDJ NT seq',
              'Insertion count', 'Insertion lengths', 'Deletion count', 'Deletion lengths',
              'VDJ cysteine count', 'CDR3 cysteine count']
 return sep.join(fields)


def _schief_output_line(seq, legacy):
    if seq is None:
        return [''] * 20
    line = []
    # line.append(seq['experiment'] if 'experiment' in seq else '')
    # line.append(seq['group'] if 'group' in seq else '')
    # line.append(seq['subject'] if 'subject' in seq else '')
    # line.append(seq['timepoint'] if 'timepoint' in seq else '')
    line.append(seq['v_gene']['gene'])
    if seq['chain'] == 'heavy':
        line.append(seq['d_gene']['gene'] if 'd_gene' in seq else '')
    line.append(seq['j_gene']['gene'])
    line.append(seq['cdr3_len'])
    line.append(seq['junc_aa'])
    line.append(seq['junc_nt'])
    line.append(100. - seq['nt_identity']['v'])
    line.append(_get_fr_identity(seq, res='nt'))
    line.append(100. - seq['aa_identity']['v'])
    line.append(_get_fr_identity(seq, res='aa'))
    line.append('yes' if 'v_ins' in seq else '')
    line.append('yes' if 'v_del' in seq else '')
    line.append(seq['vdj_aa'])
    line.append(seq['vdj_nt'])
    if 'v_ins' in seq:
        len_field = 'len' if legacy else 'length'
        line.append(len(seq['v_ins']))
        line.append('[' + ' '.join([str(i[len_field]) for i in seq['v_ins']]) + ']')
    else:
        line.append('0')
        line.append('')
    if 'v_del' in seq:
        len_field = 'len' if legacy else 'length'
        line.append(len(seq['v_del']))
        line.append('[' + ' '.join([str(i[len_field]) for i in seq['v_del']]) + ']')
    else:
        line.append('0')
        line.append('')
    if 'vdj_aa' in seq:
        line.append(seq['vdj_aa'].upper().count('C'))
    else:
        line.append('0')
        line.append('')
    if 'cdr3_aa' in seq:
        line.append(seq['cdr3_aa'].upper().count('C'))
    else:
        line.append('0')
        line.append('')
    return line

### END OF CODE FROM BRYAN

schief_csv_output(pairs, "./schief_style.csv", sep=',', legacy_abstar=False)





