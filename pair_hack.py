#!/usr/bin/env python
# filename: pair_hack.py


#
# Copyright (c) 2016 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


from abutils.core.sequence import Sequence



class Pair(object):
    '''
    Holds a pair of sequences, corresponding to HC and LC of a single mAb.

    Input is a list of dicts, with each dict containing sequence information from a single
    chain, formatted as would be returned from a query on a MongoDB database containing
    AbStar output.
    '''
    def __init__(self, seqs, name=None, h_selection_func=None, l_selection_func=None):
        self._seqs = seqs
        self._heavy = None
        self._light = None
        self._heavies = [s for s in seqs if s['chain'] == 'heavy']
        self._lights = [s for s in seqs if s['chain'] in ['kappa', 'lambda']]
        self._name = name
        self._fasta = None
        self._sample = None
        self._subject = None
        self._group = None
        self._experiment = None
        self._timepoint = None
        self._is_pair = None
        self._vrc01_like = None
        self._lineage = None
        self._select_heavy = h_selection_func
        self._select_light = l_selection_func

    def __eq__(self, other):
        return (self.heavy, self.light) == (other.heavy, other.light)

    def __ne__(self, other):
        return not self == other

    def __hash(self):
        return hash((self.heavy, self.light))


    @property
    def heavy(self):
        if self._heavy is None:
            if len(self._heavies) > 0:
                if self._select_heavy is not None:
                    self._heavy = Sequence(self._select_heavy(self._heavies))
                else:
                    self._heavy = Sequence(self._heavies[0])
            else:
                self._heavy = None
        return self._heavy

    @heavy.setter
    def heavy(self, heavy):
        self._heavy = heavy

    @property
    def light(self):
        if self._light is None:
            # self._lights = [s for s in self._seqs if s['chain'] in ['kappa', 'lambda']]
            if len(self._lights) > 0:
                if self._select_light is not None:
                    self._light = Sequence(self._select_light(self._lights))
                else:
                    self._light = Sequence(self._lights[0])
            else:
                self._light = None
        return self._light

    @light.setter
    def light(self, light):
        self._light = light

    @property
    def is_pair(self):
        if all([self.heavy is not None, self.light is not None]):
            return True
        return False

    @property
    def lineage(self):
        if self._lineage is None:
            if 'clonify' in self.heavy:
                self._lineage = self.heavy['clonify']['id']
        return self._lineage

    @property
    def vrc01_like(self):
        if self._vrc01_like is None:
            #if any([self.heavy is None, self.light is None]):
            if any([self.heavy is None]):
                self._vrc01_like = False
            else:
                self._vrc01_like = all([self.heavy['v_gene']['gene'] == 'IGHV1-2'])
                #self._vrc01_like = all([self.heavy['v_gene']['gene'] == 'IGHV1-2', self.light['cdr3_len'] == 5])
#        self._vrc01_like = True
        return self._vrc01_like


    @property
    def name(self):
        if self._name is None:
            if self.heavy is not None:
                self._name = self.heavy['seq_id']
            elif self.light is not None:
                self._name = self.light['seq_id']
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def sample(self):
        if self._sample is None:
            slist = []
            if self.experiment is not None:
                slist.append(str(self.experiment))
            if self.group is not None:
                slist.append(str(self.group))
            if self.subject is not None:
                slist.append(str(self.subject))
            if self.timepoint is not None:
                slist.append(str(self.timepoint))
            if slist:
                self._sample = '|'.join(slist)
        return self._sample

    @property
    def subject(self):
        if self._subject is None:
            if self.heavy is not None and 'subject' in self.heavy.keys():
                self._subject = self.heavy['subject']
            elif self.light is not None and 'subject' in self.light.keys():
                self._subject = self.light['subject']
        return self._subject

    @subject.setter
    def subject(self, subject):
        self._subject = subject

    @property
    def group(self):
        if self._group is None:
            if self.heavy is not None and 'group' in self.heavy.keys():
                self._group = self.heavy['group']
            elif self.light is not None and 'group' in self.light.keys():
                self._group = self.light['group']
        return self._group

    @group.setter
    def group(self, group):
        self._group = group

    @property
    def experiment(self):
        if self._experiment is None:
            if self.heavy is not None and 'experiment' in self.heavy.keys():
                self._experiment = self.heavy['experiment']
            elif self.light is not None and 'experiment' in self.light.keys():
                self._experiment = self.light['experiment']
        return self._experiment

    @experiment.setter
    def experiment(self, experiment):
        self._experiment = experiment

    @property
    def timepoint(self):
        if self._timepoint is None:
            if self.heavy is not None and 'timepoint' in self.heavy.keys():
                self._timepoint = self.heavy['timepoint']
            elif self.light is not None and 'timepoint' in self.light.keys():
                self._timepoint = self.light['timepoint']
        return self._timepoint

    @timepoint.setter
    def timepoint(self, timepoint):
        self._timepoint = timepoint


    def fasta(self, key='vdj_nt', append_chain=True):
        '''
        Returns the sequence pair as a fasta string. If the Pair object contains
        both heavy and light chain sequences, both will be returned as a single string.

        By default, the fasta string contains the 'vdj_nt' sequence for each chain. To change,
        use the <key> option to select an alternate sequence.

        By default, the chain (heavy or light) will be appended to the sequence name:

        >MySequence_heavy

        To just use the pair name (which will result in duplicate sequence names for Pair objects
        with both heavy and light chains), set <append_chain> to False.
        '''
        fastas = []
        for s, chain in [(self.heavy, 'heavy'), (self.light, 'light')]:
            if s is not None:
                c = '_{}'.format(chain) if append_chain else ''
                fastas.append('>{}{}\n{}'.format(s['seq_id'], c, s[key]))
        return '\n'.join(fastas)


def assign_pairs(seqs, name='seq_id', delim=None, delim_occurance=1, pairs_only=False):
    '''
    Assigns sequences to the appropriate mAb pair, based on the sequence name.

    Inputs:

    ::seqs:: is a list of dicts, of the format returned by querying a MongoDB containing
        Abstar output.
    ::name:: is the dict key of the field to be used to group the sequences into pairs.
        Default is 'seq_id'
    ::delim:: is an optional delimiter used to truncate the contents of the ::name:: field.
        Default is None, which results in no name truncation.
    ::delim_occurance:: is the occurance of the delimiter at which to trim. Trimming is performed
        as delim.join(name.split(delim)[:delim_occurance]), so setting delim_occurance to -1 will
        trucate after the last occurance of delim. Default is 1.
    ::pairs_only:: setting to True results in only truly paired sequences (pair.is_pair == True)
        will be returned. Default is False.

    Returns a list of Pair objects, one for each mAb pair.
    '''
    pdict = {}
    for s in seqs:
        if delim is not None:
            pname = delim.join(s[name].split(delim)[:delim_occurance])
        else:
            pname = s[name]
        if pname not in pdict:
            pdict[pname] = [s, ]
        else:
            pdict[pname].append(s)
    pairs = [Pair(pdict[n], name=n) for n in pdict.keys()]
    if pairs_only:
        pairs = [p for p in pairs if p.is_pair]
    return pairs