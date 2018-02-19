# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 10:28:37 2015

@author: zehl
"""
import odml


def fill_odml_with_xls(doc, xls, source_filename):
    """
    """
    for d in xls:
        prop_path = ':'.join([d['Path'], d['Property Name']])
        prop = doc.get_property_by_path(prop_path)
        if len(d['Value']) == 1 and len(prop.values) == 1:
            prop.value.data = d['Value'][0]
            prop.value.dtype = d['odML Data Type'][0]
            prop.value.filename = source_filename
            if len(d['Data Unit']) > 0:
                prop.value.unit = d['Data Unit'][0]
            if len(d['Data Uncertainty']) > 0:
                prop.value.uncertainty = d['Data Uncertainty'][0]

        elif len(d['Value']) > 1 and len(prop.values) == 1:
            for i, v in enumerate(d['Value']):
                if i == 0:
                    prop.value.data = v
                    prop.value.filename = source_filename
                    if len(d['Data Unit']) > 0:
                        prop.value.unit = d['Data Unit'][i]
                    if len(d['Data Uncertainty']) > 0:
                        prop.value.uncertainty = d['Data Uncertainty'][i]
                else:
                    prop.values.append(
                        odml.Value(data=v,
                                   dtype=d['odML Data Type'][i],
                                   unit=d['Data Unit'][i]
                                   if len(d['Data Unit']) > 0
                                   else '',
                                   uncertainty=d['Data Uncertainty']
                                   if len(d['Data Uncertainty']) > 0
                                   else '',
                                   filename=source_filename))

        elif len(prop.values) == len(d['Value']):
            for i, v in enumerate(d['Value']):
                prop.values[i].data = d['Value'][i]
                prop.values[i].dtype = d['odML Data Type'][i]
                prop.values[i].filename = source_filename
                if len(d['Data Unit']) > 0:
                    prop.values[i].unit = d['Data Unit'][i]
                if len(d['Data Uncertainty']) > 0:
                    prop.values[i].uncertainty = d['Data Uncertainty'][i]


def fill_odml_with_elinfo(doc, elinfo, source_filename):
    """
    """
    ff = lambda x: x.name.startswith('Electrode_')\
        and x.parent.name == 'Electrodes'
    elsecs = doc.itersections(filter_func=ff)

    for els in elsecs:
        elxinfo = elinfo[els.name]
        for elp_name in elxinfo.keys():
            els.properties[elp_name].value.data = elxinfo[elp_name]['data']
            els.properties[elp_name].value.filename = source_filename


def fill_odml_with_ba_elids(doc, ba_elids, source_filename):
    """
    """
    ff = lambda x: x.name.startswith('Electrode_')\
        and x.parent.name == 'Electrodes'
    elsecs = doc.itersections(filter_func=ff)

    for els in elsecs:
        els.properties['IDba'].value.data = ba_elids[els.name]['IDba']['data']
        els.properties['IDba'].value.filename = source_filename


def fill_odml_with_subsobj(doc, subsobj, source_filename):
    """
    """
    ff = lambda x: x.name == 'Session' and x.parent.name == 'Subsession'
    prop = [p for p in doc.iterproperties(filter_func=ff)][0]
    prop.value.data = \
        subsobj.parameters_nev['SessionName'].split('/')[-1].split('-')[0]

    ff = lambda x: x.name == 'Date' and x.parent.name == 'Subsession'
    prop = [p for p in doc.iterproperties(filter_func=ff)][0]
    prop.value.data = subsobj.parameters_nev['DateTime'].date()

    ff = lambda x: x.name == 'Weekday' and x.parent.name == 'Subsession'
    prop = [p for p in doc.iterproperties(filter_func=ff)][0]
    prop.value.data = subsobj.parameters_nev['DateTime'].strftime('%A')

    ff = lambda x: x.name == 'Time' and x.parent.name == 'Subsession'
    prop = [p for p in doc.iterproperties(filter_func=ff)][0]
    prop.value.data = \
        subsobj.parameters_nev['DateTime'].replace(microsecond=0).time()

    ff = lambda x: x.name == 'Duration' and x.parent.name == 'Subsession'
    prop = [p for p in doc.iterproperties(filter_func=ff)][0]
    prop.value.data = \
        float(subsobj.get_max_time().rescale('s').round(decimals=3))

    tasktype = doc.sections['Subsession'].properties['TaskType'].value.data
    if tasktype in ['OneCue', 'TwoCues', 'Observation']:
        ff = lambda x: x.name == 'TaskSettings' and \
            x.parent.name == 'Subsession'
        sec = [s for s in doc.itersections(filter_func=ff)][0]
        sec.properties['TotalTrialCount'].value.data = len(subsobj.trial_data)

        sec.properties['CorrectTrialCount'].value.data = \
            len(subsobj.get_trial_ids(performance_codes=[127]))

        sec.properties['GripErrorTrialCount'].value.data = \
            len(subsobj.get_trial_ids(performance_codes=[126]))

    ff = lambda x: x.name.startswith('Trial_')\
        and x.parent.name == 'TaskSettings'
    tr_secs = doc.itersections(filter_func=ff)

    trids = subsobj.trial_data[:, subsobj.trial_data_idx['trial_id']].tolist()

    for tr_sec in tr_secs:
        tr = int(tr_sec.name[-3:])
        trdata = subsobj.trial_data[trids.index(tr)]

        tr_sec.properties['ID'].value.data = tr
        tr_sec.properties['ID'].value.filename = source_filename

        pc = trdata[subsobj.trial_data_idx['performance_code']]
        tr_sec.properties['PerformanceCode'].value.data = pc
        tr_sec.properties['PerformanceCode'].value.filename = source_filename

        tt = trdata[subsobj.trial_data_idx['trial_type']]
        tr_sec.properties['TrialType'].value.data = tt
        tr_sec.properties['TrialType'].value.filename = source_filename

        for dev in subsobj.trial_events[:-2]:
            tr_sec['DigitalEvents'].properties[dev].value.data = \
                trdata[subsobj.trial_data_idx[dev]]
            tr_sec['DigitalEvents'].properties[dev].value.filename = \
                source_filename

    if len(subsobj.channel_id_nev) > 0:
        for elparams in subsobj.parameters_nev_electrodes:
            if len(elparams.keys()) > 0:
                elid = elparams['ConnectorPin'] + \
                    (elparams['ConnectorID']-1) * 32

                ff = lambda x: x.name.startswith('Electrode_%03i' % elid)\
                    and x.parent.name == 'Electrodes'
                el_secs = [s for s in doc.itersections(filter_func=ff)]
                if len(el_secs) > 0 and subsobj.is_sorted:
                    spso_sec = el_secs[0]['OfflineSpikeSorting']
                    spso_sec.properties['SpikeThreshold'].value.data = \
                        elparams['AmpThresholdLo']
                    spso_sec.properties['SpikeThreshold'].value.filename = \
                        source_filename


def fill_odml_with_aev(doc, aev, source_filename):
    """
    """
    ff = lambda x: x.name.startswith('Trial_')\
        and x.parent.name == 'TaskSettings'
    elsecs = doc.itersections(filter_func=ff)

    for els_tr in elsecs:
        tr = int(els_tr.name[-3:])
        els_tr_aev = els_tr['AnalogEvents']

        if tr in aev['trials'].keys():

            els_tr_aev.properties['UsedForceSensor'].value.data = \
                aev['sensors'][aev['trials'][tr]['trialtype']]
            els_tr_aev.properties['UsedForceSensor'].value.filename = \
                source_filename

            els_tr_aev.properties['ManuallyControlled'].value.data = \
                aev['trials'][tr]['checked']
            els_tr_aev.properties['ManuallyControlled'].value.filename = \
                source_filename

            for st in ['GripForceSignals', 'DisplacementSignal']:
                for ae in aev['trials'][tr][st].keys():
                    els_tr_aev[st].properties[ae].value.data = \
                        aev['trials'][tr][st][ae]
                    els_tr_aev[st].properties[ae].value.filename = \
                        source_filename


def fill_odml_with_alf(doc, alf, source_filename):
    """
    """
    ff = lambda x: x.name.startswith('Trial_')\
        and x.parent.name == 'TaskSettings'
    elsecs = doc.itersections(filter_func=ff)

    for els_tr in elsecs:
        els_tr.properties['ForceType'].value.data = alf[els_tr.name]
        els_tr.properties['ForceType'].value.filename = source_filename


def fill_odml_with_unitids(doc, unitids, source_filename):
    """
    """
    ff = lambda x: x.name == 'IsSpikeSorted' and x.parent.name == 'Subsession'
    prop = [p for p in doc.iterproperties(filter_func=ff)][0]
    prop.value.data = True

    ff = lambda x: x.name.startswith('Electrode_')\
        and x.parent.name == 'Electrodes'
    el_secs = doc.itersections(filter_func=ff)

    for el_sec in el_secs:
        # add sua ids
        spso_sec = el_sec['OfflineSpikeSorting']
        values = [spso_sec.properties['SUAIDs'].value.clone() for i in
                  range(len(unitids[el_sec.name]['SUAIDs']))]
        for i, v in enumerate(values):
            v.data = unitids[el_sec.name]['SUAIDs'][i]
        spso_sec.properties['SUAIDs'].values = values
        # add mua id
        spso_sec.properties['MUAID'].value.data = unitids[el_sec.name]['MUAID']
        # update noise ids
        spso_sec.properties['NoiseIDs'].values.append(
            spso_sec.properties['NoiseIDs'].value.clone())
        spso_sec.properties['NoiseIDs'].values[0].data = 0


def fill_odml_with_periods(doc, subsobj, aev):
    """
    """
    trids = sorted(subsobj.get_trial_ids())
    maxtrid = max(trids)

    for i, tr in enumerate(trids):
        ff = lambda x: x.name.startswith('Trial_%03i' % tr)\
            and x.parent.name == 'TaskSettings'
        tr_secs = [s for s in doc.itersections(filter_func=ff)]
        if len(tr_secs) > 0:
            tr_sec = tr_secs[0]
            trdata = subsobj.trial_data[trids.index(tr)]

            trev = {}
            for dev in subsobj.trial_events[:-2]:
                trev[dev] = int(trdata[subsobj.trial_data_idx[dev]] / 30.)

            for st in [['_gf', 'GripForceSignals'],
                       ['_dp', 'DisplacementSignal']]:
                if tr in aev['trials'].keys():
                    trev['OT'+st[0]] = aev['trials'][tr][st[1]]['OT']
                    trev['HS'+st[0]] = aev['trials'][tr][st[1]]['HS']
                    trev['OR'+st[0]] = aev['trials'][tr][st[1]]['OR']
                else:
                    trev['OT'+st[0]] = 0
                    trev['HS'+st[0]] = 0
                    trev['OR'+st[0]] = 0

            rt = trev['SR'] - trev['GO-ON']
            tr_sec['Periods'].properties['RT'].value.data = \
                rt if rt > 0 else -1

            if i < len(trids)-1:
                if int(trids[i+1]) == maxtrid:
                    tr_sec['Periods'].properties['ITI'].value.data = -1
                else:
                    nxt_tr = trids[i+1]
                    nxt_trdata = subsobj.trial_data[trids.index(nxt_tr)]
                    tr_sec['Periods'].properties['ITI'].value.data = \
                        int(nxt_trdata[subsobj.trial_data_idx['TS']] / 30.) - \
                        max(trev.values())

            for st in [['_gf', 'GripForceSignals'],
                       ['_dp', 'DisplacementSignal']]:
                rgt = trev['OT'+st[0]] - trev['SR']
                tr_sec['Periods'][st[1]].properties['RGT'].value.data = \
                    rgt if rgt > 0 else -1
                pt = trev['HS'+st[0]] - trev['OT'+st[0]]
                tr_sec['Periods'][st[1]].properties['PT'].value.data = \
                    pt if pt > 0 else -1
                mt = trev['HS'+st[0]] - trev['SR']
                tr_sec['Periods'][st[1]].properties['MT'].value.data = \
                    mt if mt > 0 else -1
                ht = trev['RW'] - trev['HS'+st[0]]
                tr_sec['Periods'][st[1]].properties['HT'].value.data = \
                    ht if ht > 0 else -1
                ht = trev['OR'+st[0]] - trev['RW']
                tr_sec['Periods'][st[1]].properties['ORT'].value.data = \
                    ht if ht > 0 else -1


def fill_odml_with_lfprej(doc, lfprej, subsession_name):
    """
    """
    for rejn in lfprej.keys():
        source_file = '_'.join([subsession_name, 'rej', '-'.join([
                                '%03iHz' % int(lfprej[rejn]['Lowcut']),
                                '%03iHz' % int(lfprej[rejn]['Highcut']),
                                '%02io.hdf5' % lfprej[rejn]['Order']])])
        for pn in lfprej[rejn].keys():
            if pn in ['RejElectrodes', 'RejTrials']:
                prop = doc['RejectionsLFP'][rejn].properties[pn]
                prop.value.filename = source_file
                values = [prop.value.clone() for i in
                          range(len(lfprej[rejn][pn]))]
                for i, v in enumerate(values):
                    v.data = lfprej[rejn][pn][i]
                prop.values = values
            else:
                doc['RejectionsLFP'][rejn].properties[pn].value.data = \
                    lfprej[rejn][pn]
                doc['RejectionsLFP'][rejn].properties[pn].value.filename = \
                    source_file