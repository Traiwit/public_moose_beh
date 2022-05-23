#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Tested using:
 - numpy 1.16.2
"""


import numpy as np

from bae.log_01 import msg

from topo import StaticEntities

#{ Indexable Sequence-Object
class Sequence(StaticEntities):
    """Numpy-like indexable datastructure to store sequence data/names

    The idea is to place each frame on an integer-valued time line (unitTimes),
    where one step represents a (time)unit (e.g. month, quarter, year). This
    allows a correct (time)spacing for plots.

    After __init__ you can call setupSequence to create frame/sequence-labels
    according to regularSequence.

    Example:
    ========
        >>> uTimes = [0, 12, 24] + range(25, 32)
        >>> seq = Sequence(start=('Y2017','Q4'), unit='quarterly',
        ...                 unitTimes=uTimes)
        >>> seq.setupSequence()
        >>> print seq.labels
        >>> print seq[3:-4].labels
    """

    def __init__(self, **kwargs):
        """Constructor of Sequence object.

        A Sequence object can be constructed by
         - not passing any data/keywordArguments. This results in an 'empty'
         Sequence object. Indexing this will just return an empty Seqence.
         - passing
             - 'unit', 'starts' and 'length' or
             - 'unit', 'starts' and 'unitTimes'
         keywordArguments and calling L{setupSequence} after construction.
         This will automatically setup an indexable Sequence object holding
         default sequence labels.
          >>> seq = Sequence(start=('Y2016','JUN'), unit='monthly', length=10)
          >>> seq.setupSequence()
          >>> print seq.labels[1:5]
          ['Y2016_JUL' 'Y2016_AUG' 'Y2016_SEP' 'Y2016_OCT']
         - keywordArguments having indexable values of the same (sequence)
         length wich will be assinged as indexable variables/attributes to
         Sequence object:
          >>> seq = Sequence(labelA=['a','b','c'], labelB=['1','2','3'])
          >>> print seq[[0,2]].labelA, seq[[0,2]].labelB
          ['a' 'c'] ['1' '3']

        @kwarg unitTimes: number of regular spaced sequence steps or list of
            indices of sequenceSteps with respect to given sequenceType and
            start (required)

        @kwarg unit: 'monthly', 'quarterly', 'yearly' (required)

        @kwarg starts: tuple(year, id), where year can be a string ('y2013') or
            number and id can be an index (monthly: 3-->march,
            quarterly: 2-->Q2) or string (e.g. monthly: 'mar',  quarterly:'Q2')
            (required)

        @kwarg length: number of sequence steps, required if no unitTimes given

        @kwarg unitTimes: increasing integer sequence of sequence times,
            required if no length given
        """
        super(Sequence, self).__init__()

        self.unit   = kwargs.pop('unit', None)
        self.start  = kwargs.pop('start', None)

        length      = kwargs.pop('length', None)
        unitTimes   = kwargs.pop('unitTimes', None)

        if length and unitTimes:
            msg('WARNING. Sequence definition is ambiguous. You passed ' +
                'unitTimes, so the length-argument will be overwritten.')

        if unitTimes:
            self.length = len(unitTimes)
            self.addVariable('unitTimes', unitTimes)

        if not unitTimes and length:
            self.length = length
            self.addVariable('unitTimes', range(length))

        for key, val in kwargs.iteritems():
            self.addVariable(key, val)


    def __eq__(self, other):
        '''Compares two sequence objects

        Usage:
        ======
         >>> seqA == seqB
        '''
        if not isinstance(other, self.__class__):
            return False

        for e in ['unit', 'start', 'length', 'labels', 'unitTimes']:
            if not getattr(self, e) == getattr(other, e):
                return False
        else:
            return True

    def setupSequence(self):
        ''' Name labels using regularSequence function
        '''
        if self.unitTimes is None and self.length is not None:
            #only length is given
            self.addVariables('unitTimes', range(self.length))

        # check if everything is available
        if not any([d is None
                    for d in [self.unit, self.start, self.unitTimes]]):

            labels = regularSequence(self.unitTimes, self.unit, self.start)
            self.addVariable('labels', labels)
        else:
            raise ValueError('The attributes unit, start and'
                             ' unitTimes/length are needed to setup sequence.')
#} #End Sequence



###############################################################################
#{Frames

### Some usefull functions for default framenaming
def defaultNamesFromIds(ids):
    '''Creates defaultNames from given ids: [(2,4),] --> ['F2004',]
    '''
    return np.array(['F%d%03d' % tuple(id) for id in ids])


def defaultIdsFromNames(names, defaultStep=2, frameDigits=3):
        '''Parses frameNames to ids
        @param defaultStep: if names are type of [F001, F002, ...] this stepId
            will be used (--> [(1,1), (1,2), ]). If step is included in
            names (e.g. ['F2001', 'F2002',...] ) the stepId will be parsed.
            default = 2
        @param frameDigits: specifies number of digits used for frame
            default = 3
        '''

        import re
        def parser(s):
            nums = re.findall(r'\d+', s)
            if len(nums) == 1:
                if int(nums[0]) > 10**frameDigits - 1:
                    # Type of F2014
                    return (int(nums[0][0]), int(nums[0][1:]))
                else:
                    # Type of F014
                    return (defaultStep, int(nums[0]))
            elif len(nums) == 2:
                #Type of Step2Frame014
                return (int(nums[0]), int(nums[1]))
            else:
                raise ValueError("Can't parse a valid frameId-tuple from "
                                 "'%s'" % s)

        frameIds = map(parser, names)

        return frameIds


###maybe: Frames --> TimePoints? (TimeFrames?)
class Frames(StaticEntities):
    """array of time points for timedependent data.

    Example:
    ========
        >>> frameIds = [(2,d) for d in range(10)]
        >>> ff = Frames(ids=frameIds) # will create default (frame)names
        >>> customLabels = ['custom_%d'%ii for ii in range(10)]
        >>> ff.addTimeVariable('custom', customLabels)
        >>> ff.regularSequence(unit='quarterly', start=(2018,'Q2'))
        >>> ffSub = ff[3:5]
        >>> print ffSub.ids
        >>> print ffSub.names
        >>> print ffSub.seq.labels
        >>> print ffSub.custom

    """

    def __init__(self, **kwargs):
        """
        Frames can be initialized by keywordArguments arguments. Common entrys
        are:
         - 'ids': a list of (step,frame)-tuples: [(2,0), (2,1), ...],
         (tuples have to be unique and in increasing order)
         - 'names': a list of frameNames ['F2000','F2001', ...]
         - 'seq': an index-able L{Sequence}-object
        and/or as .

        Example:
         >>> ids   = [(2,1),(2,2)]
         >>> names = ['F2001','F2002']
         >>> seq   = Sequence(unit='quarterly', start(2019,'Q1'), length=2)
         >>> fr1 = Frames(ids = ids, names = names, seq = seq)
         >>> #or
         >>> fr2 = Frames()
         >>> fr2.addVariable(ids=ids)
         >>> fr2.addVariable(names=names)
         >>> fr2.addVariable(seq=seq)
         >>> # if names or ids not given, __init__ will try to autocomplete
         >>> # these variables
         >>> fr4 = Frames(ids=ids)
         >>> fr5 = Frames(names=names)
         >>> # fr4 and fr5 will hold the same data as fr1,... because names
         >>> # are given in a standard form here
        """
        StaticEntities.__init__(self, **kwargs)

        if ('ids' in self.variables) and not ('names' in self.variables):
            self._defaultNamesFromIds()

        if ('names' in self.variables) and not ('ids' in self.variables):
            try:
                self._defaultIdsFromNames()

            except ValueError as e:
                msg("WARNING! While initializing Frame object " +
                    "the FrameIds (not passed here) couldn't be "+
                    "determined only from frameNames: \n %s" % str(e))


    def idsAsTupleList(self):
        """Returns the frameIds as a list of tuples, e.g.: [(2,1), (2,4), ...]

        @note: requires the attributes 'ids' to exist
        """
        return map(tuple, self.ids)


    def addVariable(self, name, values):
        """Adds a mutable attribute to self, e.g. a topo field (i.e. point
        list) or a time series.

        If name=='ids' then check if ids are unique and in ascending order.

        @param values: must be a numpy array or convertible to such. Or
        it must have a __getitem__ and a __len__ function. __getitem__ must
        then accept arguments like a numpy-array.

        @Note: Silently overwrites an already existing attribute of the same
           name!
        """
        if name == 'ids':
            # 1st check if ascending are valid and unique
            values = np.asarray(values)
            sVals, idx = np.unique(values, axis=0, return_index=True)

            if (not idx.shape[0] == values.shape[0] or
                not np.all(idx==np.arange(idx.shape[0]))):
                raise ValueError('Sequence-ids must be unique and of '+
                                 'ascending order')
        StaticEntities.addVariable(self, name, values)


    def idxFromIds(self, ids, strict=True):
        '''Returns the index of a given sequence of frameIds. Both, input
        ids and stored ids have to be unique, otherwise an error is raised.

        @param ids: sequence of frameids (e.g. [(2,0),...]) or numpy array of
            shape (N,2)

        @param strict: if True, not all requeset ids has to be found
            (default: True)

        @note: requires the attributes 'ids' to exist
        '''
        ids = np.asarray(ids)

        uIds, fIdx, cnt = np.unique(np.vstack((self.ids, ids)),
                                    return_index=True,
                                    return_counts=True,
                                    axis=0)

        cntNotFound = uIds.shape[0] - self.ids.shape[0]

        if cntNotFound:
            text = '%d requested ids can not be found' % cntNotFound
            if strict:
                raise ValueError(text)
            else:
                msg('WARNING! ' + text)

        if np.any(cnt>2):
            raise ValueError('self.ids and/or input-ids have not been unique')

        return fIdx[(cnt == 2) & (fIdx < self.ids.shape[0])]


    def idxFromNames(self, names):
        '''Returns indices for given sequence of frameNames

        @note: requires the attributes 'names' to exist
        '''
        if isinstance(names, str):
            names = [names,]

        try:
            index = [self.names.index(name) for name in names]
        except ValueError:
            raise ValueError("Can't find frame with id %s." % str(id))

        return index


    def regularSequence(self, **kwargs):
        '''Creates a regular indexable L{Sequence} object using
        the L{regularSequence} function and appends it to indexable attributes
        of Frames.

        Example
        =======
         >>> f = Frames([(2,ii) for ii in range(10)])
         >>> # initialy 3yearly (12quarter) steps, than quarterly
         >>> uTimes = [0, 12, 24] + range(25, 32) #length of 10
         >>> f.regularSequence(start=('Y2017','Q4'),
         ...                  unit='quarterly',
         ...                  unitTimes=uTimes)
         >>> print f[:-1].seq.labels

        @kwarg unit: 'monthly', 'quarterly', 'yearly'

        @kwarg start: tuple(year, id), where year can be a string ('y2013') or
            a number and id can be an index (monthly: 3-->march,
            quarterly: 2-->Q2) or string (e.g. monthly: 'mar',  quarterly:'Q2')

        @kwarg length: number of sequence steps, required if no unitTimes given

        @kwarg unitTimes: increasing integer sequence of sequence times,
            required if no length given
        '''
        if 'unitTimes' not in kwargs.keys():
            kwargs['length'] = len(self)
        else:
            if not len(kwargs['unitTimes']) == len(self):
                raise ValueError(('Length of unitTimes (%d) must match'
                                  ' the number of Frames (%d)') %
                                 (len(kwargs['unitTimes']), len(self)))

        seq = Sequence(**kwargs)

        if not hasattr(self, 'labels'):
            seq.setupSequence()

        self.addVariable('seq', seq)


    def _defaultNamesFromIds(self):
        '''Creates defaultNames from given self.ids: [(2,4),] --> ['F2004',]

        @note: requires the attributes 'ids' to exist
        '''
        names = defaultNamesFromIds(self.ids)

        self.addVariable('names', names)


    def _defaultIdsFromNames(self, defaultStep=2, frameDigits=3):
        '''Parses frameNames to ids

        @param defaultStep: if names are type of [F001, F002, ...] this stepId
            will be used (--> [(1,1), (1,2), ]). If step is included in
            names (e.g. ['F2001', 'F2002',...] ) the stepId will be parsed.

        @param frameDigits: specifies number of digits used for frame

        @note: requires the attributes 'names' to exist
        '''
        ids = defaultIdsFromNames(self.names, defaultStep=defaultStep,
                                  frameDigits=frameDigits)

        self.addVariable('ids', ids)
#} #End Frames

###############################################################################
#{Sequence

### sequence generation
def regularSequence(unitTimes, unit, starts):
    '''Creates a regular sequence

    Example:
    ========
     >>> # create a 'quarterly' sequence, where the last steps are not
     >>> # equaly spaced
     >>> unitPoints = range(15) + range(15, 22, 2) + [25, 31]
     >>> seq = regularSequence(unitPoints, 'quarterly', (2019, 'Q2'))

    @param unitTimes: number of regular spaced sequence steps or list of
        indices of sequenceSteps with respect to given sequenceType and start

    @param unit: 'monthly', 'quarterly', 'yearly'

    @param starts: tuple(year, id), where year can be a string ('y2013') or
        number and id can be an index (monthly: 3-->march,
        quarterly: 2-->Q2) or string (e.g. monthly: 'mar',  quarterly:'Q2')
    '''

    if isinstance(unitTimes, (int,long)):
        unitTimes = range(unitTimes)
    else:
        if not all([b>a for a,b in zip(unitTimes[:-1], unitTimes[1:])]):
            raise ValueError("Sequence indices must have ascending order")

    if unit.lower() in 'monthly':
        lookup = ["JAN","FEB","MAR","APR","MAY","JUN",
                  "JUL","AUG","SEP","OCT","NOV","DEC"]

    elif unit.lower() in 'quarterly':
        lookup = ["Q%d" % (d+1) for d in range(4)]

    elif unit.lower() in 'yearly':
        lookup = ['']

    try:
        startYear, sIndex = starts

    except TypeError as e:
        if unit.lower() in 'yearly':
            sIndex = 0

        else:
            print 'Error-origin', str(e)
            raise ValueError('starts argument must be a ' +
                             '(year, index/string)-tuple ' +
                             'for unit %s' % unit)

    if isinstance(sIndex, str):
        try:
            sIndex = [l.lower() for l in lookup].index(sIndex.lower())
        except ValueError:
            raise ValueError("StartIndex for sequence can't be parsed." +
                             " The passed identifier %s is not in lookup %s"
                             % (sIndex, str(lookup)))

    else:
        sIndex -= 1

    try:
        lookup[sIndex]
    except IndexError:
        raise IndexError('The index %s is not valid for unit %s' %
                         (str(sIndex + 1), unit))

    if isinstance(startYear, str):
        startYear = int( startYear.lower().replace('y','') )

    nSeq = []
    for kk in unitTimes:
        offYear, idx = divmod( (kk + sIndex), len(lookup) )
        year = startYear + offYear
        part = lookup[idx]
        if part:
            nSeq.append('Y%s_%s' % (year, part))
        else:
            nSeq.append(year)

    nSeq = np.array(nSeq)

    return nSeq




###############################################################################

### Testfuctions
def _test_Frames():

    # test frames
    ee = Frames()

    mm = Frames(ids=[(2,d) for d in range(10)])

    ff = Frames(ids=[(2,d) for d in range(10)],names=
                 ['F%03d' % d for d in range(10)])



    ll = ff.idsAsTupleList()[2:5]
    ii = ff.idxFromIds(ll)


    # test Sequence
    uTimes = [0, 12, 24] + range(25, 32)
    seq = Sequence(start=('Y2017','Q4'), unit='quarterly', unitTimes=uTimes)
    seq.setupSequence()
    msg( seq.labels )
    msg( seq[3:-4].labels )

    # test Sequence implementation in Frames
    ff.regularSequence(unit='quarterly', start=(2018,'Q2'))



if __name__ == '__main__':
    #_test_Frames()
    print 'No syntax errors.'
