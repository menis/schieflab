#!/usr/bin/python3
from abutils.core import sequence
from abutils.core import pair
from abutils.utils.alignment import global_alignment, muscle
from pair_hack import Pair as Pair_hack
import json
import sys, os
import plots
from datetime import date

### Additional for plotting
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
from abtools import color
#from abtools.alignment import global_alignment, muscle
import plots

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import Align

from collections import OrderedDict
###

# Globals
VH12 = 'IGHV1-2'
VK320 = 'IGKV3-20'
VK133 = 'IGKV1-33'
VL214 = 'IGLV2-14'
VK15 = 'IGKV1-5'

def munge_description(x):
    tokens = x.split(" ")
    num_of_tokens = len(tokens)
    if num_of_tokens != 3:
        return dict(cluster_fraction=tokens[0])
    else:
        return dict(cluster_fraction=tokens[1], cluster_confidence=tokens[2])

seqs = []
## Parse args
if '-file' in sys.argv:
    TargetFile = sys.argv[sys.argv.index('-file') + 1]
    basename = os.path.splitext(os.path.basename(TargetFile))[0]
    with open(TargetFile) as f:
        for line in f:
            if line.strip() == "":
                continue
            d = json.loads(line.strip())
            # Bryan has a pair object
            seq = sequence.Sequence(d)
            seqs.append(seq)

force_all_heavy_as_vrc01class = False
blankrun = False

if '-forcevrc01' in sys.argv:
    force_all_heavy_as_vrc01class = True

if '-blankrun' in sys.argv:
    blankrun = True

colortouse = '#45bc70'
if '-plotcolor' in sys.argv:
    color = sys.argv[sys.argv.index('-plotcolor') + 1]
    if color == "orange":
        colortouse = '#f99248'
    elif color == "green":
        colortouse = '#45bc70'
    else:
        colortouse = str(color).strip()

expanded = False
tag = ""
if '-expanded' in sys.argv:
    expanded = True
    tag = "_expanded"

species = 'human'
if '-species' in sys.argv:
    species = sys.argv[sys.argv.index('-species') + 1]

# Residue masks to include in VRC01-class mutation counting
vh12mask = False
vk320mask = False
vk133mask = False
vl214mask = False
vk15mask = False

masks = OrderedDict()
if '-mask' in sys.argv:
    maskfile = sys.argv[sys.argv.index('-mask') + 1]
    with open(maskfile) as f:
        startfound = False
        currentchain = ""
        currentdict = OrderedDict()
        position = 1
        for line in f:
            l = line.strip()
            if l == "":
                continue
            elif l.split()[0] == "chain":
                if startfound:
                    # found a new chain, save the current dict into the masks
                    masks[currentchain] = currentdict
                    currentdict = OrderedDict()
                    currentchain = l.split()[1]
                    position = 1
                else:
                    # found a new chain but this is the first one, reset the dict just in case
                    # and note the chain
                    currentdict = OrderedDict()
                    currentchain = l.split()[1]
                    startfound = True
            else:
                # all lines should be int he form of "Q yes" or "Q no". exit if anything else
                lineparts = l.split()
                shoudlidie = False
                errorstring = ""

                if len(lineparts) != 2:
                    shoudlidie = True
                    errorstring += " More than two items per line;"
                else:
                    if len(lineparts[0]) > 1:
                        shoudlidie = True
                        errorstring += " First part should be a single character aa code;"
                    if lineparts[1] not in ("yes", "no"):
                        shoudlidie = True
                        errorstring += " Second part should be a yes or no;"

                if shoudlidie:
                    exit("ERROR in the mask file. Line: " + l + errorstring)
                else:
                    currentdict[position] = lineparts[1]
                    position += 1
        # If there is still something in the dict at the end that means we need to store the last chain
        # before moving on
        if len(currentdict) != 0:
            masks[currentchain] = currentdict
            currentdict = OrderedDict()

        # Set the flags
        if VH12 in masks:
            vh12mask = True
        if VK320 in masks:
            vk320mask = True
        if VK133 in masks:
            vk133mask = True
        if VL214 in masks:
            vl214mask = True
        if VK15 in masks:
            vk15mask = True

# munge the ids and create a dictionary from which we construct
# the cellid and the chainAnnotation.
for s in seqs:
    s.dictionary = munge_description(str(s['description']))
    s.fraction = s.dictionary['cluster_fraction'] if 'cluster_fraction' in s.dictionary else 'Unknown'
    s.confidence = s.dictionary['cluster_confidence'] if 'cluster_confidence' in s.dictionary else 'Unknown'

# order the sequences by their munged id.
# seqordered = [seq for seq in sorted(seqs, key = lambda x: munge(x.id), reverse = True)]
seqordered = [seq for seq in sorted(seqs, key=lambda x: x.id, reverse=True)]

# Next we want to 'cluster' sequences with the same cellid
cells = dict.fromkeys([s.id for s in seqordered], None)

# dict.fromkeys ends up with each nested dict as a reference of the same thing, so we loop instead and explicitly create different nested dicts.
for k in cells.keys():
    cells[k] = {'heavy': [], 'light': []}

for s in seqordered:
    if (s['chain'] == 'heavy'):
        cells[s.id]['heavy'].append(s)
    else:
        cells[s.id]['light'].append(s)

# Pair sorted sequences based on munged id.
# Unpaired sequences aren't exported. Should they be?
pairs = []
unpaired = []
# Within a 'cluster' of sequences we want to pair each heavy and each light chain sequence.
for cell in cells:
    if (len(cells[cell]['heavy']) == 0) and (len(cells[cell]['light']) > 0):
        unpaired.extend([pair.Pair([s]) for s in cells[cell]['light']])
    if (len(cells[cell]['light']) == 0) and (len(cells[cell]['heavy']) > 0):
        unpaired.extend([pair.Pair([s]) for s in cells[cell]['heavy']])
    if (len(cells[cell]['light']) > 0) and (len(cells[cell]['heavy']) > 0):
        # for each heavy pair it with each light.
        for heavy in cells[cell]['heavy']:
            for light in cells[cell]['light']:
                pairs.append(pair.Pair([heavy, light]))

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
    shared = []
    total = []
    mutationpositions = []
    make_dir(os.path.dirname(output_file))
    header = _get_schief_output_header(sep)
    output = [header, ]
    for p in sorted(pairs, key=lambda x: _get_name(x)):
        name = _get_name(p)
        line = [name, ]
        line += _get_pair_metadata(p)
        vrc01class = False
        if p.vrc01_like:
            vrc01class = True
        elif p.heavy is not None:
            hp = Pair_hack([p.heavy])
            if hp.vrc01_like and force_all_heavy_as_vrc01class:
                vrc01class = True
        line += _schief_output_line(p.heavy, legacy_abstar, vrc01class, shared, total, mutationpositions)
        line += _schief_output_line(p.light, legacy_abstar, vrc01class, ch="light")
        output.append(sep.join([str(l) for l in line]))
    open(output_file, 'w').write('\n'.join(output))
    return shared, total, mutationpositions


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
    fields = ['Sequence ID', 'Experiment', 'Group', 'Subject', 'Timepoint', 'Heavy Cluster Fraction', 'Heavy Cluster Confidence', 'VH gene', 'VH Alternate:score',
              'DH gene', 'DH Alternate:score', 'JH gene', 'JH Alternate:score',
              'H_CDR3 length', 'H_CDR3 nt', 'H_CDR3 aa',
              'H_Junction AA', 'H_Junction NT seq', 'H_V NT muts', 'H_V AA muts', 'H_D NT muts', 'H_D AA muts', 'H_J NT muts', 'H_J AA muts',  '% VH mutation (NT)', '% H FR mutation (NT)',
              '% VH mutation (AA)', '% H FR mutation (AA)', 'H_VRC01-like Mutations', 'VH insertions', 'VH deletions',
              'H_VDJ AA seq', 'H_VDJ NT seq', 'H_Insertion count', 'H_Insertion lengths:position',
              'H_Deletion count', 'H_Deletion lengths:position', 'H_VDJ cysteine count', 'H_CDR3 cysteine count', 'Light Cluster Fraction', 'Light Cluster Confidence',
              'VL gene', 'VL Alternate:score', 'JL gene', 'JL Alternate:score', 'L_CDR3 length', 'L_CDR3 nt', 'L_CDR3 aa',
              'Junction AA', 'Junction NT seq', 'V NT muts', 'V AA muts', 'J NT muts', 'J AA muts',
              '% VL mutation (NT)', '% L FR mutation (NT)', '% VL mutation (AA)', '% L FR mutation (AA)', 'L_VRC01-like Mutations',
              'VL insertions', 'VL deletions', 'L_VDJ AA seq', 'L_VDJ NT seq',
              'L_Insertion count', 'L_Insertion lengths:position', 'L_Deletion count', 'L_Deletion lengths:position',
              'L_VDJ cysteine count', 'L_CDR3 cysteine count']
    return sep.join(fields)


def _schief_output_line(seq, legacy, pairisvrc01class=False, s=None, t=None, m=None, ch="heavy"):
    if seq is None:
        if ch == "heavy":
            return [''] * 34
        else:
            return [''] * 30
    line = []
    line.append(seq.fraction)
    line.append(seq.confidence)
    line.append(seq['v_gene']['gene'])
    line.append(_get_alternates(seq['v_gene']['others'], j_gene=False))
    if seq['chain'] == 'heavy':
        line.append(seq['d_gene']['gene'] if 'd_gene' in seq else '')
        line.append(_get_alternates(seq['d_gene']['others'], j_gene=False) if 'd_gene' in seq else '')
    line.append(seq['j_gene']['gene'])
    line.append(_get_alternates(seq['j_gene']['others'], j_gene=True))
    line.append(seq['cdr3_len'] if 'cdr3_len' in seq else '')
    line.append(seq['cdr3_nt'].upper() if 'cdr3_nt' in seq else '')
    line.append(seq['cdr3_aa'].upper() if 'cdr3_aa' in seq else '')
    line.append(seq['junc_aa'])
    line.append(seq['junc_nt'])
    line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in seq['var_muts_nt']['muts']) if 'var_muts_nt' in seq else '')
    line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in seq['var_muts_aa']['muts']) if 'var_muts_aa' in seq else '')
    if seq['chain'] == 'heavy':
        line.append('Data Not Available')
        line.append('Data Not Available')
        # D Mutations are currently not avaialbe in Abstar
        # line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in seq['d_muts_nt']['muts'])  if 'd_muts_nt' in seq else '')
        # line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in seq['d_muts_aa']['muts'])  if 'd_muts_aa' in seq else '')
    line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in seq['join_muts_nt']['muts']) if 'join_muts_nt' in seq else '')
    line.append('|'.join(str(i['was']) + str(i['position']) + str(i['is']) for i in seq['join_muts_aa']['muts']) if 'join_muts_aa' in seq else '')
    line.append(100. - seq['nt_identity']['v'])
    line.append(_get_fr_identity(seq, res='nt'))
    line.append(100. - seq['aa_identity']['v'])
    line.append(_get_fr_identity(seq, res='aa'))

    if pairisvrc01class:
        justvgene = seq['v_gene']['aa_sequence']
        name = seq['seq_id']
        trimmed = sequence.Sequence(justvgene)
        trimmed['seq_id'] = name
        trimmed['vdj_aa'] = justvgene
        if seq['chain'] == 'heavy':
            vrc01_class, total = vrc01_class_mutation_count([trimmed], vgene_only=True)
            muts = vrc01_class_mutation_positions([trimmed], vgene_only=True)
            s.extend(vrc01_class)
            t.extend(total)
            m.extend(muts)
            line.append(vrc01_class[0])
        else:
            vgeneforthischain = pickclosestvgene(trimmed, species, seq['v_gene']['gene'])
            vrc01_class, total = vrc01_class_mutation_count_light_chain([trimmed], vgene_only=True, vgene=vgeneforthischain)
            #muts = vrc01_class_mutation_positions([trimmed], vgene_only=True)
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
    line.append(seq['vdj_aa'].upper().count('C') if 'vdj_aa' in seq else '')
    line.append(seq['cdr3_aa'].upper().count('C') if 'cdr3_aa' in seq else '')
    return line


### END OF CODE FROM BRYAN

### Additional code from Bryan to count VRC01 Mutations
###

def vrc01_class_mutation_count(seqs, vgene_only=True):
    input_seqs = [sequence.Sequence([s['seq_id'], s['vdj_aa']]) for s in seqs]

    # get VRC01-class sequences
    if(not expanded):
        vrc01_seqs = get_vrc01_class_sequences(vgene_only=vgene_only)
    else:
        vrc01_seqs = get_expanded_vrc01_class_sequences(vgene_only=vgene_only)

    # get glVRC01 sequence
    glvrc01 = get_vrc01_germline_sequence(vgene_only=vgene_only)

    return identifyvrc01muts(input_seqs, vrc01_seqs, glvrc01, VH12)

def get_mask_string(germline, maskforthischain):
    position = 1
    mstring = ""
    for c in germline:
        if c == '-':
            mstring += "0"
        else:
            if maskforthischain[position] == "yes":
                mstring += "1"
            else:
                mstring += "0"
            position += 1

    return mstring


def identifyvrc01muts(input_seqs, vrc01_seqs, glvrc01, genecall):
    import re
    regex = re.compile("[a-zA-Z]")

    shared = []
    total = []
    glvrc01_name = glvrc01.id
    vrc01_names = [s.id for s in vrc01_seqs]

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
            aln_seq.seq = aln_seq.seq[index[0]:(index[-1]+1)]
            aln_gl.seq = aln_gl.seq[index[0]:(index[-1]+1)]
            for p in aln_vrc01s:
                p.seq = p.seq[index[0]:(index[-1]+1)]

        # If there is a VH12 Mask. Generate a mask string 1 for include and 0 otherwise
        chainmaskstring = "x" * len(aln_gl.seq)
        usemask = False
        if genecall == VH12 and vh12mask:
            chainmaskstring = get_mask_string(aln_gl.seq, masks[VH12])
            usemask = True
        elif genecall == VK320 and vk320mask:
            chainmaskstring = get_mask_string(aln_gl.seq, masks[VK320])
            usemask = True
        elif genecall == VK133 and vk133mask:
            chainmaskstring = get_mask_string(aln_gl.seq, masks[VK133])
            usemask = True
        elif genecall == VL214 and vl214mask:
            chainmaskstring = get_mask_string(aln_gl.seq, masks[VL214])
            usemask = True

        # Count mutations
        _total = sum([_s != g for _s, g in zip(str(aln_seq.seq), str(aln_gl.seq)) if g != '-'])

        total.append(_total)
        all_shared = {}
        for vrc01 in aln_vrc01s:
            _shared = []
            for q, g, v, cmask in zip(str(aln_seq.seq), str(aln_gl.seq), str(vrc01.seq), str(chainmaskstring)):
                if g == '-' and v == '-':
                    _shared.append(False)
                elif q == v and q != g:
                    if not usemask:
                        _shared.append(True)
                    else:
                        if cmask == "1":
                            _shared.append(True)
                        else:
                            _shared.append(False)
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

def vrc01_class_mutation_count_light_chain(seqs, vgene_only=True, vgene=VK320):
    input_seqs = [sequence.Sequence([s['seq_id'], s['vdj_aa']]) for s in seqs]

    shared = []
    total = []

    # get VRC01-class sequences
    if(not expanded):
        vrc01_seqs = get_vrc01_class_sequences(chain=vgene, vgene_only=vgene_only)
    else:
        vrc01_seqs = get_expanded_vrc01_class_sequences(vgene_only=vgene_only)

    # get glVRC01 sequence
    glvrc01 = get_vrc01class_germline_lights(vgene)

    return identifyvrc01muts(input_seqs, vrc01_seqs, glvrc01, vgene)

def pickclosestvgene(seq, species='human', assignment=VK320):
    if species == 'human' and assignment in (VK320, VK133, VL214, VK15):
        return assignment
    else:
        ## Need to find the closest match
        lightchains_gl = get_vrc01class_germline_lights('all')
        currentscore = -1000
        closestlight = VK320
        for chain in lightchains_gl:
            aligner = Align.PairwiseAligner()
            aligner.match_score = 3
            aligner.mismatch_score = -2
            aligner.open_gap_score = -22
            aligner.extend_gap_score = -1
            score = aligner.score(seq.sequence, chain.sequence)
            if score > currentscore:
                currentscore = score
                closestlight = chain.id

        return closestlight

def vrc01_class_mutation_positions(seqs, vgene_only=True):
    data = []
    input_seqs = [sequence.Sequence([s['seq_id'], s['vdj_aa']]) for s in seqs]
    input_names = [s.id for s in input_seqs]
    # get VRC01-class sequences
    if(not expanded):
        hiv_seqs = get_vrc01_class_sequences()
    else:
        hiv_seqs = get_expanded_vrc01_class_sequences()

    all_hiv_names = [s.id for s in hiv_seqs]
    # MSA
    seqs_for_alignment = input_seqs + hiv_seqs
    seqs_for_alignment.append(get_vrc01_germline_sequence(vgene_only=vgene_only))
    aln = muscle(seqs_for_alignment)
    aln_seqs = [seq for seq in aln if seq.id in input_names]
    aln_gl = [seq for seq in aln if seq.id == 'glVRC01'][0]
    aln_mins = [seq for seq in aln if seq.id in ['minVRC01', 'min12A21']]
    aln_hiv = [seq for seq in aln if seq.id in all_hiv_names]
    for seq in aln_seqs:
        seq_data = []
        for i, (s, g) in enumerate(zip(str(seq.seq), str(aln_gl.seq))):
            # if g == '-' and s == '-':
            if g == '-':
                continue
            min_residues = [seq[i] for seq in aln_mins]
            vrc01_residues = [seq[i] for seq in aln_hiv]
            if s == '-':
                seq_data.append(0)
            elif s == g:
                seq_data.append(0)
            elif s != g and s in min_residues:
                seq_data.append(2)
            elif s != g and s in vrc01_residues:
                seq_data.append(3)
            elif s != g and s not in vrc01_residues:
                seq_data.append(1)
            else:
                seq_data.append(0)
        data.append(np.asarray(seq_data))
    return np.asarray(data)

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
        vk3_20_light = [('VRC01',
                         'EIVLTQSPGTLSLSPGETAIISCRTSQYGSLAWYQQRPGQAPRLVIYSGSTRAAGIPDRFSGSRWGPDYNLTISNLESGDFGVYYCQQYEFFGQGTKVQVDIKR'),
                        ('PGV04',
                         'EIVLTQSPGTLSLSPGETASLSCTAASYGHMTWYQKKPGQPPKLLIFATSKRASGIPDRFSGSQFGKQYTLTITRMEPEDFARYYCQQLEFFGQGTRLEIRR'),
                        ('3BNC60',
                         'DIQMTQSPSSLSARVGDTVTITCQANGYLNWYQQRRGKAPKLLIYDGSKLERGVPARFSGRRWGQEYNLTINNLQPEDVATYFCQVYEFIVPGTRLDLKRTVAA')]
        vk1_33_light = [('VRC-CH31',
                         'DIQMTQSPSSLSASLGDRVTITCQASRGIGKDLNWYQQKAGKAPKLLVSDASTLEGGVPSRFSGSGFHQNFSLTISSLQAEDVATYFCQQYETFGQGTKVDIK'),
                        ('12A12',
                         'DIQMTQSPSSLSASVGDRVTITCQAGQGIGSSLQWYQQKPGKAPKLLVHGASNLHRGVPSRFSGSGFHTTFSLTISGLQRDDFATYFCAVLEFFGPGTKVEIKRTVAAPSV')]
        vl2_14_light = [('PGV20', 'QSALTQPPSVSGSPGQSITLSCTGASTSVAWYQQYADKAPRLIVFDGNKRPSDISSRFSGSQSGGTASLTISGLQSEDEAYYHCNAFEFFGGGTKLTVL')]

        vk1_5_light = [('N60P23', 'CFVQSQSPSTLYASVGDKITITCRSSQSGWKAWYQQKPGKAPKLLIYKSSILETGVPSRFIGSDSGTEFTLTISSLQPDDFATYYCQHFETFG')]

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
    if chain == 'heavy':
        seqs = heavy
    elif chain == VK320:
        seqs = vk3_20_light
    elif chain == VK133:
        seqs = vk1_33_light
    elif chain == VL214:
        seqs = vl2_14_light
    elif chain == VK15:
        seqs = vk1_5_light
    else:
        seqs = light

    if only_include is not None:
        if type(only_include) in [str, unicode]:
            only_include = [only_include, ]
        seqs = [s for s in seqs if s[0] in only_include]
    return [sequence.Sequence(s) for s in seqs]

def get_expanded_vrc01_class_sequences(chain='heavy', vgene_only=True, only_include=None):
    heavy = [('VRC01', 'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTR'),
             ('PGV04', 'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCAR'),
             ('VRC-CH31', 'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCAR'),
             ('3BNC60', 'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCAR'),
             ('12A12', 'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCAR'),
             ('PGV20', 'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCAR'),
             ('PCIN63-71Ia','QVQLVQSGAEVKKPGASVRVSCKASGYTFNSCLIHWWRQAPGQGLQWMAWINPLHGAVNYAHQFQGRITVTRDTSIDTAYMELRGLRSDDTATYYCTR')]
    light = []

    seqs = heavy if chain == 'heavy' else light
    if only_include is not None:
        if type(only_include) in [str, unicode]:
            only_include = [only_include, ]
        seqs = [s for s in seqs if s[0] in only_include]
    return [sequence.Sequence(s) for s in seqs]

def _get_mutations(seqs, standard):
    mutations = []
    for seq in seqs:
        aln = global_alignment(seq, target=standard,
            matrix='blosum62', gap_open=-15, gap_extend=-1)
        mutations.extend(_parse_mutations(aln))
    return mutations

def get_vrc01_class_mutations():

    if(not expanded):
        vrc01_class = [s.sequence for s in get_vrc01_class_sequences()]
    else:
        vrc01_class = [s.sequence for s in get_expanded_vrc01_class_sequences()]
    glvrc01 = get_vrc01_germline_sequence().sequence
    return list(set(_get_mutations(vrc01_class, glvrc01)))

def get_expanded_vrc01_class_mutations():
    vrc01_class = [s.sequence for s in get_expanded_vrc01_class_sequences()]
    glvrc01 = get_vrc01_germline_sequence().sequence
    return list(set(_get_mutations(vrc01_class, glvrc01)))

def get_vrc01_germline_sequence(vgene_only=True):
    if vgene_only:
        gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR')

        # Hack to get to work with JH sequences. Issue is that alignment is stuid sometimes at the N-term
        #gl_vrc01 = ('glVRC01', 'PGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR')
    else:
        gl_vrc01 = ('glVRC01',
                    'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGKNSDYNWDFQHWGQGTLVTVSS')
    return sequence.Sequence(gl_vrc01)

def get_vrc01class_germline_lights(clone=VK320):

    vk3_20_gl = ('IGKV3-20',
                  'EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSP')
    vk1_33_gl = ('IGKV1-33',
                 'DIQMTQSPSSLSASVGDRVTITCQASQDISNYLNWYQQKPGKAPKLLIYDASNLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLP')
    vl2_14_gl = ('IGLV2-14',
                  'QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYEVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTL')

    vk1_5_gl = ('IGKV1-5', 'DIQMTQSPSTLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCQQYNSYS')

    all = [('IGKV3-20',
                  'EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSP'),
           ('IGKV1-33',
                 'DIQMTQSPSSLSASVGDRVTITCQASQDISNYLNWYQQKPGKAPKLLIYDASNLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLP'),
           ('IGLV2-14',
            'QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYEVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTL')
           ]

    if clone == VK320:
        lights = sequence.Sequence(vk3_20_gl)
    elif clone == VK133:
        lights = sequence.Sequence(vk1_33_gl)
    elif clone == VL214:
        lights = sequence.Sequence(vl2_14_gl)
    elif clone == VK15:
        lights = sequence.Sequence(vk1_5_gl)
    else:
        lights = [sequence.Sequence(s) for s in all]

    return lights

def _parse_mutations(aln):
    muts = []
    tpos = 0
    for q, t in zip(aln.aligned_query, aln.aligned_target):
        if t == '-':
            continue
        tpos += 1
        if any([q == '-', q == t]):
            continue
        mut = '{}{}'.format(tpos, q)
        muts.append(mut)
    return muts

#Copied from https://github.com/chris-rands/CrUtils/blob/master/scripts/faTranslateBioPython.py
#https://github.com/chris-rands
def pad_seq(seq):
    """Pad sequence to multiple of 3 with Ns"""
    return {0: seq, 1: seq+'NN', 2: seq+'N'}[len(seq) % 3]


## Generate trendline


vrc01_class_mutations = get_vrc01_class_mutations()
# Donor 5684
d5684_mutations = []
with open('/Users/menis/src/schieflabscripts/schieflab/synthetic_antibodies/d5684_mutations.txt') as f:
    for line in f:
        d5684_mutations.append(line.strip().split())
# there are actually 1920 synthetic sequences from donor 5684,
# but we only need the first 1000
d5684_mutations = d5684_mutations[:1000]


# Donor 6471
d6471_mutations = []
with open('/Users/menis/src/schieflabscripts/schieflab/synthetic_antibodies/d6471_mutations.txt') as f:
    for line in f:
        d6471_mutations.append(line.strip().split())

# combine synthetic mutations from both donors
#synth_mutations = d5684_mutations + d6471_mutations

mean = [0, ]
ci_min = [0, ]
ci_max = [0, ]

# iterate temporally through the mutations - take the first mutation,
# then the first two, etc...
for i in range(1, 21):
    vrc01_class = []

    #need to split this out to know which mutations come from where
    for sm in d5684_mutations:
        # only sample sequences that have at least i mutations
        if len(sm) < i:
            continue
        mutsfound = sum([m in vrc01_class_mutations for m in sm[:i]])
        vrc01_class.append(mutsfound)
        if(blankrun):
            print("d5684", i, mutsfound)

    for sm in d6471_mutations:
        # only sample sequences that have at least i mutations
        if len(sm) < i:
            continue
        mutsfound = sum([m in vrc01_class_mutations for m in sm[:i]])
        vrc01_class.append(mutsfound)
        if (blankrun):
            print("d6471", i, mutsfound)


    # calculate the mean and 95% CI for the frequency of VRC01-class mutations

    n, min_max, _mean, var, skew, kurt = stats.describe(vrc01_class)
    std = np.sqrt(var)
    _ci_min, _ci_max = stats.norm.interval(0.95, loc=_mean, scale=std / np.sqrt(len(vrc01_class)))
    mean.append(_mean)
    ci_min.append(_ci_min)
    ci_max.append(_ci_max)

if(blankrun):
    #don't need to continue with the rest
    sys.exit()

# in order to center the squares in our 2D histogram, we need to offset
# each datapoint by 0.5
rand_xs = [x + 0.5 for x in range(len(mean))]
mean = [m + 0.5 for m in mean]
ci_min = [c + 0.5 for c in ci_min]
ci_max = [c + 0.5 for c in ci_max]

total_vrc01_mut = []
total_muts = []
mut_positions = []

vr, tot, mu = schief_csv_output(pairs, "./"+basename+"_pairs_"+date.today().strftime("%d%B%Y")+tag+".csv", sep=',', legacy_abstar=False)
total_vrc01_mut.extend(vr)
total_muts.extend(tot)
mut_positions.extend(mu)

# Note that if force_all_heavy_as_vrc01class is not enabled the vr and tot counts will be zero for unpaired chains
vr, tot, mu = schief_csv_output(unpaired, "./"+basename+"_unpaired_"+date.today().strftime("%d%B%Y")+tag+".csv", sep=",", legacy_abstar=False)
total_vrc01_mut.extend(vr)
total_muts.extend(tot)
mut_positions.extend(mu)

## Print graphs

# minval=0.05 ensures that the lightest boxes are still distinguishable from background
cmap = color.truncate_colormap(color.cmap_from_color(colortouse), minval=0.05)

# plot the frequency of VRC01-class mutation in the GT8 immunization group
f, ax = plt.subplots()
if(len(total_muts) > 0 and len(total_vrc01_mut) > 0):
 plots.shared_mutation_2dhist(total_muts, total_vrc01_mut, cmap, ax, show_values=False)

# plot the frequency of random VRC01-class mutation, with 95% CIs
plots.fill_between_steps(ax, rand_xs, ci_max, ci_min)
plt.step(rand_xs, mean, where='mid', alpha=0.9, c='k', linewidth=0.5)

# Extra bells and whistles
ax.plot([0, 19], [0, 19], transform=ax.transAxes, linestyle='--', color='red', linewidth=1)
points = [[0, 1], [0, 15], [15, 16]]
t = plt.Polygon(points, color='#f2f2f2')
plt.gca().add_patch(t)


plt.tight_layout()
plt.savefig("./"+basename+"_vrc01countplot_"+date.today().strftime("%d%B%Y")+tag+".pdf")
plt.savefig("./"+basename+"_vrc01countplot_"+date.today().strftime("%d%B%Y")+tag+".svg")


schief_csv_output(pairs+unpaired, "./"+basename+"_paired_and_unpaired_"+date.today().strftime("%d%B%Y")+tag+".csv", sep=",", legacy_abstar=False)

cmap = ListedColormap(['#F5F5F5', '#080808', '#0DABE6', '#0DABE6'])
if(len(mut_positions) > 0):
    plots.pixel_plot(np.asarray(mut_positions),
                 cmap, figfile="./"+basename+"_mutationdistribution_"+date.today().strftime("%d%B%Y")+tag+".pdf"
                )

