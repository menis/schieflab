#!/opt/conda/bin/python
from abutils.core import sequence
from abutils.core import pair
from abutils.utils.alignment import global_alignment, muscle
from pair_hack import Pair as Pair_hack
import json
import sys, os

seqs = []
with open(sys.argv[1]) as f:
    for line in f:
        d = json.loads(line.strip())
        seq = sequence.Sequence(d)
        seqs.append(seq)

# order the sequences by their munged id.
# seqordered = [seq for seq in sorted(seqs, key = lambda x: munge(x.id), reverse = True)]
seqordered = [seq for seq in sorted(seqs, key=lambda x: x.id, reverse=True)]

# Pair sorted sequences based on munged id.
# Unpaired sequences aren't exported. Should they be?
pairs = []
unpaired = []

first = None
nxt = None

while seqordered:
   if(first is None):
       first = seqordered.pop()
   if(seqordered):
       if(nxt is None):
           nxt = seqordered.pop()
       if munge(first.id) == munge(nxt.id):
           pairs.append(pair.Pair([first,nxt]))
           first = None
           nxt = None
       else:
           unpaired.append(pair.Pair([first]))
           first = nxt
           nxt = None
   else:
       unpaired.append(pair.Pair([first]))
       first = None

print("Processed", len(pairs), "pairs of sequences")
print("and ", len(unpaired), " unpaired sequences")


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
        # line += _get_pair_metadata(p)
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
    experiment = seq['experiment'] if 'experiment' in seq else ''
    group = seq['group'] if 'group' in seq else ''
    subject = seq['subject'] if 'subject' in seq else subject
    timepoint = seq['timepoint'] if 'timepoint' in seq else subject
    subject = seq.dictionary['ptid'] if 'ptid' in seq.dictionary else subject
    timepoint = seq.dictionary['visit'] if 'visit' in seq.dictionary else timepoint
    return [experiment, group, subject, timepoint]


def _get_fr_identity(seq, res='nt'):
    len_field = 'region_len_nt' if res == 'nt' else 'region_len_aa'
    mut_field = 'region_muts_nt' if res == 'nt' else 'region_muts_aa'
    regions = ['fr1', 'fr2', 'fr3']
    length = sum([seq[len_field][region] for region in regions])
    muts = sum([len(seq[mut_field][region]['muts']) for region in regions])
    return 100. * muts / length


def _get_alternates(seq, j_gene=False):
    alternates = []
    for entry in seq:
        if j_gene:
            alternates.append(str(entry['full']) + ":" + str(entry['score']))
        else:
            alternates.append(str(entry['full']) + ":" + str(entry['assigner_score']))
    if len(alternates) > 0:
        return '|'.join(map(str, alternates))
    else:
        return ' '


def _get_schief_output_header(sep):
    fields = ['Sequence ID', 'Experiment', 'Group', 'Subject', 'Timepoint', 'VH gene', 'VH Alternate:score',
              'DH gene', 'DH Alternate:score', 'JH gene', 'JH Alternate:score',
              'CDR3 length', 'CDR3 nt', 'CDR3 aa',
              'Junction AA', 'Junction NT seq', 'NT muts', 'AA muts', '% VH mutation (NT)', '% FR mutation (NT)',
              '% VH mutation (AA)', '% FR mutation (AA)', 'VRC01-like Mutations', 'VH insertions', 'VH deletions',
              'VDJ AA seq', 'VDJ NT seq', 'Insertion count', 'Insertion lengths:position',
              'Deletion count', 'Deletion lengths:position', 'VDJ cysteine count', 'CDR3 cysteine count',
              'VL gene', 'VL Alternate:score', 'JL gene', 'JL Alternate:score', 'CDR3 length', 'CDR3 nt', 'CDR3 aa',
              'Junction AA', 'Junction NT seq', 'NT muts', 'AA muts',
              '% VL mutation (NT)', '% FR mutation (NT)', '% VL mutation (AA)', '% FR mutation (AA)',
              'VL insertions', 'VL deletions', 'VDJ AA seq', 'VDJ NT seq',
              'Insertion count', 'Insertion lengths:position', 'Deletion count', 'Deletion lengths:position',
              'VDJ cysteine count', 'CDR3 cysteine count']
    return sep.join(fields)


def _schief_output_line(seq, legacy):
    if seq is None:
        return [''] * 27
    line = []
    line.append(seq['experiment'] if 'experiment' in seq else '')
    line.append(seq['group'] if 'group' in seq else '')
    line.append(seq.dictionary['ptid'] if 'ptid' in seq.dictionary else '')
    line.append(seq.dictionary['visit'] if 'visit' in seq.dictionary else '')
    #    line.append(seq['subject'] if 'subject' in seq else '')
    #    line.append(seq['timepoint'] if 'timepoint' in seq else '')
    line.append(seq['v_gene']['gene'])
    line.append(_get_alternates(seq['v_gene']['others'], j_gene=False))
    if seq['chain'] == 'heavy':
        line.append(seq['d_gene']['gene'] if 'd_gene' in seq else '')
        line.append(_get_alternates(seq['d_gene']['others'], j_gene=False) if 'd_gene' in seq else '')
    line.append(seq['j_gene']['gene'])
    line.append(_get_alternates(seq['j_gene']['others'], j_gene=True))
    line.append(seq['cdr3_len'])
    line.append(seq['cdr3_nt'].upper())
    line.append(seq['cdr3_aa'].upper())
    line.append(seq['junc_aa'])
    line.append(seq['junc_nt'])
    line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in
                         seq['var_muts_nt']['muts'] + seq['join_muts_nt']['muts']))
    line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in
                         seq['var_muts_aa']['muts'] + seq['join_muts_aa']['muts']))
    line.append(100. - seq['nt_identity']['v'])
    line.append(_get_fr_identity(seq, res='nt'))
    line.append(100. - seq['aa_identity']['v'])
    line.append(_get_fr_identity(seq, res='aa'))
    if seq['chain'] == 'heavy':
        hacky_pair = Pair_hack([seq])
        if hacky_pair.vrc01_like:
            vrc01_class, total = vrc01_class_mutation_count([hacky_pair.heavy], vgene_only=True)
            line.append(vrc01_class[0])
        else:
            line.append('')

    line.append('yes' if 'v_ins' in seq else '')
    line.append('yes' if 'v_del' in seq else '')
    line.append(seq['vdj_aa'])
    line.append(seq['vdj_nt'])
    if 'v_ins' in seq:
        len_field = 'len' if legacy else 'length'
        line.append(len(seq['v_ins']))
        line.append('[' + ' '.join([str(i[len_field]) + ":" + str(i['position']) for i in seq['v_ins']]) + ']')
    else:
        line.append('0')
        line.append('')
    if 'v_del' in seq:
        len_field = 'len' if legacy else 'length'
        line.append(len(seq['v_del']))
        line.append('[' + ' '.join([str(i[len_field]) + ":" + str(i['position']) for i in seq['v_del']]) + ']')
    else:
        line.append('0')
        line.append('')
    line.append(seq['vdj_aa'].upper().count('C'))
    line.append(seq['cdr3_aa'].upper().count('C'))
    return line


### END OF CODE FROM BRYAN

### Additional code from Bryan to count VRC01 Mutations
###

def vrc01_class_mutation_count(seqs, vgene_only=True):
    input_seqs = [sequence.Sequence([s['seq_id'], s['vdj_aa']]) for s in seqs]

    shared = []
    total = []

    # get VRC01-class sequences
    vrc01_seqs = get_vrc01_class_sequences(vgene_only=vgene_only)
    vrc01_names = [s.id for s in vrc01_seqs]

    # get glVRC01 sequence
    glvrc01 = get_vrc01_germline_sequence(vgene_only=vgene_only)
    glvrc01_name = glvrc01.id

    import re
    regex = re.compile("[a-zA-Z]")

    # identify VRC01-class mutations
    for s in input_seqs:
        alignment_seqs = [s] + vrc01_seqs + [glvrc01]
        aln = muscle(alignment_seqs)

        aln_seq = [seq for seq in aln if seq.id == s.id][0]
        aln_gl = [seq for seq in aln if seq.id == glvrc01_name][0]
        aln_vrc01s = [seq for seq in aln if seq.id in vrc01_names]

        # Strip '-' off the front of strings based on input string
        # match = re.search(regex,str(aln_seq.seq))
        index = [m.start() for m in re.finditer(regex, str(aln_seq.seq))]

        # Logic if no match is found...
        if (len(index) > 1):
            aln_seq.seq = aln_seq.seq[index[0]:index[len(index) - 1]]
            aln_gl.seq = aln_gl.seq[index[0]:index[len(index) - 1]]
            for p in aln_vrc01s:
                p.seq = p.seq[index[0]:index[len(index) - 1]]

        # Count mutations
        _total = sum([_s != g for _s, g in zip(str(aln_seq.seq), str(aln_gl.seq)) if g != '-'])

        total.append(_total)
        all_shared = {}
        for vrc01 in aln_vrc01s:
            _shared = []
            for q, g, v in zip(str(aln_seq.seq), str(aln_gl.seq), str(vrc01.seq)):
                if g == '-' and v == '-':
                    _shared.append(False)
                elif q == v and q != g:
                    _shared.append(True)
                else:
                    _shared.append(False)
            all_shared[vrc01.id] = _shared
        any_shared = 0
        for pos in zip(*all_shared.values()):
            if any(pos):
                any_shared += 1
        shared.append(any_shared)
        # print("Seq: %20s, Total: %2d, Shared: %2d" % (str(aln_seq.id), _total, any_shared))
        # print("Seq: "+str(aln_seq.id)+" Total: "+str(_total)+", Shared: "+str(any_shared))
    return shared, total


def get_vrc01_class_sequences(chain='heavy', vgene_only=True, only_include=None):
    if vgene_only:
        heavy = [('VRC01',
                  'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTR'),
                 ('PGV04',
                  'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCAR'),
                 ('VRC-CH31',
                  'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCAR'),
                 ('3BNC60',
                  'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCAR'),
                 ('12A12',
                  'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCAR'),
                 ('PGV20',
                  'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCAR')]
        light = []
    else:
        heavy = [('VRC01',
                  'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS'),
                 ('PGV04',
                  'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCARQKFYTGGQGWYFDLWGRGTLIVVSS'),
                 ('VRC-CH31',
                  'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCARAQKRGRSEWAYAHWGQGTPVVVSS'),
                 ('3BNC60',
                  'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCARQRSDFWDFDVWGSGTQVTVSS'),
                 ('12A12',
                  'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCARDGSGDDTSWHLDPWGQGTLVIVSA'),
                 ('PGV20',
                  'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCARRMRSQDREWDFQHWGQGTRIIVSS')]
        light = []
    seqs = heavy if chain == 'heavy' else light
    if only_include is not None:
        if type(only_include) in [str, unicode]:
            only_include = [only_include, ]
        seqs = [s for s in seqs if s[0] in only_include]
    return [sequence.Sequence(s) for s in seqs]


def get_vrc01_germline_sequence(vgene_only=True):
    if vgene_only:
        #        gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR')

        # Hack to get to work with JH sequences. Issue is that alignment is stuid sometimes at the N-term
        gl_vrc01 = ('glVRC01', 'PGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR')
    else:
        gl_vrc01 = ('glVRC01',
                    'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGKNSDYNWDFQHWGQGTLVTVSS')
    return sequence.Sequence(gl_vrc01)


schief_csv_output(pairs, "./schief_style_pairs.csv", sep=',', legacy_abstar=False)
schief_csv_output(unpaired, "./schief_style_unpaired.csv", sep=",", legacy_abstar=False)