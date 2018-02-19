# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 12:13:32 2015

@author: zehl
"""
import odml
import datetime
import itertools
import numpy as np
import quantities as pq

def print_summary(odmldoc):
    """
    Prints a summarized version of the loaded odml-document.

    Args:
        --

    Printout:
        Author: Author
        LastChanges: Date
        Section_Name (Section_Type):
            Property_Name:
                data: data values
                dtype: data value type
                unit: data value unit
            ...
            Subsection_Name (Subsection_Type):
                Property_Name:
                    data: data values
                    dtype: data value type
                    unit: data value unit
                ...

    Notes:
        * For values of dtype text only the first fife words are printed.
        * For values of dtype list only the first and last two values are
          printed.
    """
    # FIXME not working yet
    print 'Author:', odmldoc.author
    print 'LastChanges:', odmldoc.date

    spaths = np.sort([s.get_path() for s in odmldoc.itersections()])
    for sp in spaths:
        sobj = odmldoc.find_by_path(sp)
        sind = '\t' * sp.count('/')

        print ' '.join([sind, sp.split('/')[-1],'('+sobj.type+'):'])

        for p in sobj.properties:
            data, dtype, unit = self.mdata.get_values(p)

            if len(data) > 4:
                dshort = [data[i] for i in [0, 1, -2, -1]]
                dshort.insert(2, '...')
                data = ', '.join([str(s) for s in dshort])
            else:
                data = ', '.join([str(s) for s in data])

            if dtype in ['text', 'string']:
                if len(data.split(' ')) > 6:
                    dshort = data.split(' ')[0:6]
                    dshort.append('...')
                data = ' '.join([str(s) for s in dshort])

            print ' '.join([sind+'\t', p.name+':'])
            print ' '.join([sind+'\t\t', 'data:', data])
            print ' '.join([sind+'\t\t', 'dtype:', str(dtype)])
            print ' '.join([sind+'\t\t', 'unit:', str(unit)])


def default_value(dtype, unit=None):
    """
    Returns an odML Value with default data of the given DType

    Args:
        dtype (odml.DType):
            DType of the default data
        unit (string, optional):
            Unit of the default data

    Returns:
        (odml.value.BaseValue):
            odML Value with default data, dtype and unit

    Note:
        Valid DTypes: odml.DType.text, odml.DType.string, odml.DType.person,
        odml.DType.url, odml.DType.binary, odml.DType.boolean, odml.DType.int,
        odml.DType.date, odml.DType.datetime
    """
    if dtype in [odml.DType.text, odml.DType.string, odml.DType.person,
                 odml.DType.binary]:
        data = '-'
    elif dtype == odml.DType.url:
        data = 'file://-'
    elif dtype == odml.DType.boolean:
        data = False
    elif dtype == odml.DType.int:
        data = -1
    elif dtype == odml.DType.float:
        data = -1.0
    elif dtype == odml.DType.date:
        data = datetime.date(1111, 11, 11)
    elif dtype == odml.DType.datetime:
        data = datetime.datetime(1111, 11, 11, 11, 11, 11)
    elif dtype == odml.DType.time:
        data = datetime.time(11, 11, 11)
    else:
        raise TypeError('Unknown DType')

    return odml.Value(data=data, dtype=dtype, unit=unit)


def save_odml(doc, save_to):
    """
    Saves the odML Document as odml file to specified location

    Args:
        doc (odml.doc.BaseDocument):
            odML Document to save
        save_to (string):
            Absolute path to the desired file location
    """
    odml.tools.xmlparser.XMLWriter(doc).write_file(save_to)


def load_odml(load_from):
    """
    Loads odML Document from  specified location

    Args:
        load_from (string):
            Absolute path to the desired file location

    Returns:
        (odml.doc.BaseDocument)
    """
    doc = odml.tools.xmlparser.load(load_from)

    return doc


def get_TaskType(doc):
    """
    Returns the task type

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (str):
            task type
    """
    sec = doc['Project']['TaskDesigns']
    output = sec.properties['UsedDesign'].value.data

    return output


def get_TrialCount(doc, trialtype=None, performance_code=None):
    """
    Returns a list of trials ids

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialtype (str or int):
            If stated, returns only count of trials with given trial type
        performance_code (int):
            If stated, returns only count of trials with given performance code

    Returns:
        (int):
            Number of specified trials
    """
    sec = doc['Subsession']['TaskSettings']

    if performance_code == 127:
        output = sec.properties['CorrectTrialCount'].value.data

    elif performance_code == 126:
        output = sec.properties['GripErrorTrialCount'].value.data

    elif performance_code in [64, 96, 112, 120, 124]:
        subsec = sec['TrialTypeSettings']

    else:
        output = sec.properties['TotalTrialCount'].value.data

    # TODO: extend to trial types and other performance codes

    return output


def get_TrialIDs(doc):
    """
    Returns a list of trials ids

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (list of int):
            Trial id list
    """
    output = []

    ff = lambda x: x.name == 'TaskSettings'
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name.startswith('Trial_')
    for trsec in sec.itersections(filter_func=ff):
        FF = lambda x: x.name == 'ID'
        output.append(
            [p for p in trsec.iterproperties(filter_func=FF)][0].value.data)

    return sorted(output)


def get_PerformanceCode(doc, trialid, code=True):
    """
    Returns the performance of the monkey in the given trial either as code or
    abbreviation.

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial
        code (boolean):
            If True (default), integer code of trial performance is returned
            If False, abbreviation of trial performance is returned

    Returns:
        (int or string):
            performance code or abbreviation for wanted trial
    """
    ff = lambda x: x.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'PerformanceCode'
    output = [p for p in sec.iterproperties(filter_func=ff)][0].value.data

    if code:
        return output

    else:
        ff = lambda x: x.name == 'PerformanceCodes'
        sec = [s for s in doc.itersections(filter_func=ff)][0]

        ff = lambda x: x.name == 'pc_%i' % output
        output = [p for p in sec.iterproperties(filter_func=ff)][0].value.data

        return output


def get_TrialType(doc, trialid, code=True):
    """
    Returns trial type (code or abbreviation) for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial
        code (boolean):
            If True (default), integer code of trial type is returned
            If False, string abbreviation of trial type is returned

    Returns:
        (int or str):
            trial type for wanted trial
    """
    ff = lambda x: x.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'TrialType'
    output = [p for p in sec.iterproperties(filter_func=ff)][0].value.data

    if code:
        return output

    else:
        ff = lambda x: x.name == 'TrialTypeCodes'
        sec = [s for s in doc.itersections(filter_func=ff)][0]

        ff = lambda x: x.value.data == output
        output = [p for p in sec.iterproperties(filter_func=ff)][0].name

        return output


def get_OccurringTrialTypes(doc, code=True):
    """
    Returns all occurring trial types (code or abbreviations)

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        code (boolean):
            If True, integer code of trial type is returned
            If False, string abbreviation of trial type is returned

    Returns:
        (list of int or str):
            list of occurring trial types
    """
    trial_id_list = get_TrialIDs(doc)

    output = np.unique([get_TrialType(doc, trid, code=code) for trid in
                        trial_id_list]).tolist()

    return output


def get_ForceType(doc, trialid):
    """
    Returns force type (LF or HF) for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (str):
            force type for wanted trial
    """
    ff = lambda x: x.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'ForceType'
    output = [p for p in sec.iterproperties(filter_func=ff)][0].value.data

    return output


def get_DigitalEvent(doc, event, trialid):
    """
    Returns time of specified digital event for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        event (string):
            Abbreviation of wanted digital trial event (possible)
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            Time of wanted digital event for wanted trial
    """
    ff = lambda x: x.name == 'DigitalEvents' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == event
    output = [p for p in sec.iterproperties(filter_func=ff)][0].value.data
    unit = [p for p in sec.iterproperties(filter_func=ff)][0].value.unit

    return pq.Quantity(output, pq.CompoundUnit(unit))


def get_DigitalEvent_abbr(doc):
    """
    Returns a list of abbreviations of all digital events.

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (list of strings):
            all digital event abbreviations
    """
    ff = lambda x: x.name == 'TrialSettings' and \
        x.parent.name == 'Project'
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'DigitalEvents'
    prop = [p for p in sec.iterproperties(filter_func=ff)][0]

    output = [v.data for v in prop]

    return output


def get_AnalogEvent(doc, event, trialid):
    """
    Returns time of specified analog event for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        event (string):
            Abbreviation of wanted analog trial event (possible)
        trialid (int):
            ID of wanted trial

    Returns:
        (dict):
            keys (str):
                States which analog signal was used to calculate the wanted
                analog event ('DisplacementSignal' or 'GripForceSignals')
            values (quantities.quantity.Quantity):
                Time of wanted analog event for wanted trial for
                corresponding analog signal
    """
    output = {}

    ff = lambda x: x.name == 'AnalogEvents' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == event
    props = [p for p in sec.iterproperties(filter_func=ff)]

    for p in props:
        output[p.parent.name] = pq.Quantity(p.value.data,
                                            pq.CompoundUnit(p.value.unit))

    return output


def get_AnalogEvent_abbr(doc):
    """
    Returns a list of abbreviations of all analog events.

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (list of strings):
            all analog event abbreviations
    """
    ff = lambda x: x.name == 'TrialSettings' and \
        x.parent.name == 'Project'
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'AnalogEvents'
    prop = [p for p in sec.iterproperties(filter_func=ff)][0]

    output = [v.data for v in prop]

    return output


def get_RT(doc, trialid):
    """
    Returns reaction time for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            reaction time for wanted trial
    """
    ff = lambda x: x.name == 'Periods' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'RT'
    output = [p for p in sec.iterproperties(filter_func=ff)][0].value.data
    unit = [p for p in sec.iterproperties(filter_func=ff)][0].value.unit

    return pq.Quantity(output, pq.CompoundUnit(unit))


def get_ITI(doc, trialid):
    """
    Returns consecutive inter-trial-interval for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            inter-trial-interval for wanted trial
    """
    ff = lambda x: x.name == 'Periods' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'ITT'
    output = [p for p in sec.iterproperties(filter_func=ff)][0].value.data
    unit = [p for p in sec.iterproperties(filter_func=ff)][0].value.unit

    return pq.Quantity(output, pq.CompoundUnit(unit))


def get_RGT(doc, trialid):
    """
    Returns consecutive reach-grasp time for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            reach-grasp time for wanted trial
    """
    output = {}

    ff = lambda x: x.name == 'Periods' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'RGT'
    props = [p for p in sec.iterproperties(filter_func=ff)]

    for p in props:
        output[p.parent.name] = pq.Quantity(p.value.data,
                                            pq.CompoundUnit(p.value.unit))

    return output


def get_PT(doc, trialid):
    """
    Returns consecutive pull time for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            pull time for wanted trial
    """
    output = {}

    ff = lambda x: x.name == 'Periods' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'PT'
    props = [p for p in sec.iterproperties(filter_func=ff)]

    for p in props:
        output[p.parent.name] = pq.Quantity(p.value.data,
                                            pq.CompoundUnit(p.value.unit))

    return output


def get_MT(doc, trialid):
    """
    Returns consecutive movement time for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            movement time for wanted trial
    """
    output = {}

    ff = lambda x: x.name == 'Periods' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'MT'
    props = [p for p in sec.iterproperties(filter_func=ff)]

    for p in props:
        output[p.parent.name] = pq.Quantity(p.value.data,
                                            pq.CompoundUnit(p.value.unit))

    return output


def get_HT(doc, trialid):
    """
    Returns consecutive hold time for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            hold time for wanted trial
    """
    output = {}

    ff = lambda x: x.name == 'Periods' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'HT'
    props = [p for p in sec.iterproperties(filter_func=ff)]

    for p in props:
        output[p.parent.name] = pq.Quantity(p.value.data,
                                            pq.CompoundUnit(p.value.unit))

    return output


def get_ORT(doc, trialid):
    """
    Returns consecutive object-release time for wanted trial

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialid (int):
            ID of wanted trial

    Returns:
        (quantities.quantity.Quantity):
            object-release time for wanted trial
    """
    output = {}

    ff = lambda x: x.name == 'Periods' and \
        x.parent.name == 'Trial_%03i' % trialid
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'ORT'
    props = [p for p in sec.iterproperties(filter_func=ff)]

    for p in props:
        output[p.parent.name] = pq.Quantity(p.value.data,
                                            pq.CompoundUnit(p.value.unit))

    return output


def get_cnd(doc):
    """
    Returns condition code

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (int):
            condition code
    """
    ff = lambda x: x.name == 'TrialTypeSettings'
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'Condition'
    output = [p for p in sec.iterproperties(max_depth=0,
                                            filter_func=ff)][0].value.data

    return output


def get_trialids_trty(doc, trialtype):
    """
    Returns a list of trials ids which have the given trial type

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialtype (int or str):
            trial type of wanted trials

    Returns:
        (list of int):
            Trial id list with the given trial type
    """
    trialids = get_TrialIDs(doc)

    if isinstance(trialtype, int):
        code = True
    else:
        code = False

    output = []
    for trid in trialids:
        if get_TrialType(doc, trid, code) == trialtype:
            output.append(trid)

    return output


def get_trialids_pc(doc, performance_code):
    """
    Returns a list of trials ids which have the given performance code

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project
        trialtype (int or str):
            trial type of wanted trials

    Returns:
        (list of int):
            Trial id list with the given trial type
    """
    trialids = get_TrialIDs(doc)

    if isinstance(performance_code, int):
        code = True
    else:
        code = False

    output = []
    for trid in trialids:
        if get_PerformanceCode(doc, trid, code) == performance_code:
            output.append(trid)

    return output


def get_trialids_trtyseq(doc):
    """
    Returns a dictionary containing lists with trial ids which are part of a
    trial sequence with same trial type which follows on a trial having a
    different trial type (trial type transition) or on an error trial (error
    transition)

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (dict of list of lists of int):
            keys (tuple):
                (tt1, tt2)
    """
    # get trial types
    trialids = get_TrialIDs(doc)

    trialtypes = dict(
        (trid, get_TrialType(doc, trid, code=False)) for trid in trialids)

    octrty = get_OccurringTrialTypes(doc, code=False)

    output = {}
    for tto in itertools.product(octrty, octrty):
        if tto[0] != tto[1] and tto[1] > 0:
            output[tto] = []

    for idx, trid in enumerate(trialids):
        tr_tt_now = trialtypes[trid]

        if idx < len(trialids)-2:
            tr_tt_next = trialtypes[trialids[idx+1]]

        if tr_tt_now != tr_tt_next and tr_tt_next > 0:
            tto = (tr_tt_now, tr_tt_next)
            i = 1
            tt_in_seq = []
            if idx+i < len(trialids)-2:
                while trialtypes[trialids[idx+i]] == tr_tt_next:
                    tt_in_seq.append(trialids[idx+i])
                    i = i+1

            if len(tt_in_seq) > 0:
                if len(output[tto]) == 0:
                    output[tto] = [[] for t in range(len(tt_in_seq))]
                elif len(output[tto]) < len(tt_in_seq):
                    extend_by = len(tt_in_seq)-len(output[tto])
                    output[tto].extend([[] for t in range(extend_by)])

                for s, tr in enumerate(tt_in_seq):
                    output[tto][s].append(tr)

    return output


def get_TrialType_order(doc):
    """
    Returns the order of trial types for a subsession

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (dict):
            keys (string): 'OrderGrip' and 'OrderForce'
            values (string): 'random' or 'block'

    Note:
        'block' order of trial types means that either grip and/or force types
        are repeated in block of trials of certain size.
        'random' order of trial types means that either grip and/or force types
        are changing according to a random pattern generator (TODO: check
        details)
    """
    sec = doc['Subsession']['TaskSettings']['TrialTypeSettings']
    output = {p: sec.properties[p].value.data for p in ['OrderGrip',
                                                        'OrderForce']}

    return output


def get_TrialType_randomness(doc):
    """
    """
    trno_ctr = get_TrialCount(doc, performance_code=127)
    trno_gertr = get_TrialCount(doc, performance_code=126)

    # get trial types
    trialids = get_TrialIDs(doc)

    trialtypes = dict(
        (trid, get_TrialType(doc, trid, code=False)) for trid in trialids)

    octrty = get_OccurringTrialTypes(doc, code=False)
    if 'n.d.' in octrty:
        octrty.remove('n.d.')

    trtyseq = {}
    for tto in itertools.product(octrty, octrty):
        trtyseq[tto] = []

    for idx, trid in enumerate(trialids):
        tr_tt_now = trialtypes[trid]

        if idx < len(trialids)-2:
            tr_tt_next = trialtypes[trialids[idx+1]]
            tto = (tr_tt_now, tr_tt_next)
            if tto in trtyseq.keys():
                trtyseq[tto].append(trid)

    # get number of trials which have a trial type and are followed by a trial
    # which has no trial type
    trno_fby_tt0 = 0
    for k in trtyseq.keys():
        if k[0] > 0 and k[1] == 0:
            trno_fby_tt0 += len(trtyseq[k])

    # number of detectable trial types not followed by trial type 0
    trno_tt_inseq = trno_ctr + trno_gertr - trno_fby_tt0

    # expected number of consecutive (real) trial type orders
    expnum4randdist = float(trno_tt_inseq) / (len(octrty)**2)

    # create matrix out of trial type order results
    ttomatrix = np.zeros((len(octrty), len(octrty)))
#    test = list(np.zeros((len(octt_excltt0), len(octt_excltt0)), dtype=str))

    # create index dictionary for trial type order matrix
    ttomatrix_idx = dict([(tt, i) for (i, tt) in enumerate(octrty)])

    # fill trial type order matrix and test list
    for tto in trtyseq.keys():
        if 0 not in tto:
            ttomatrix[ttomatrix_idx[tto[0]], ttomatrix_idx[tto[1]]] = \
                np.round(len(trtyseq[tto]) / expnum4randdist, decimals=1)
#            test[ttomatrix_idx[tto[0]]][ttomatrix_idx[tto[1]]] = '%i/%i' % tto
#    test = np.asarray(test)  # transform test list to array


def is_standard(doc):
    """
    Returns if subsession was performed under standard settings

    Args:
        doc (odml.doc.BaseDocument):
            odML Document of reach-to-grasp project

    Returns:
        (boolean):
            True, is standard settings were used
    """
    ff = lambda x: x.name == 'TaskSettings'
    sec = [s for s in doc.itersections(filter_func=ff)][0]

    ff = lambda x: x.name == 'StandardSettings'
    output = [p for p in sec.iterproperties(max_depth=0,
                                            filter_func=ff)][0].value.data

    return output