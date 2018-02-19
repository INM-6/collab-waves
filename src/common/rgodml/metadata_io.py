# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 11:56:19 2015

@author: zehl
"""
import os
import csv
import glob
import xlwt
import xlrd
import types
import scipy
import datetime

import numpy as np

import rgodml.blackrock_tree as blackrock
import rgodml.project_tree as project
import rgodml.setup_tree as setup
import rgodml.subsession_tree as subsession
import rgodml.subject_tree as subject
import rgodml.preprocessing_tree as preprocessing

#import rgodml.odml_utils
import rgodml.odml_io

import h5py_wrapper.wrapper as h5py


def odml_default_value():
    """
    Returns an odML Value with default data of the given DType

    Args:
        dtype (odml.DType):
            DType of the default data

    Returns:
        (odml.value.BaseValue):
            odML Value with default data, dtype and unit

    Note:
        Valid DTypes: odml.DType.text, odml.DType.string, odml.DType.person,
        odml.DType.url, odml.DType.binary, odml.DType.boolean, odml.DType.int,
        odml.DType.date, odml.DType.datetime
    """
    defaults = ['-', 'file://-', 'False', -1, -1.0,
                datetime.date(1111, 11, 11).isoformat(),
                datetime.datetime(1111, 11, 11, 11, 11, 11).isoformat(),
                datetime.time(11, 11, 11).isoformat()]

    return defaults


def write_odml2dict(doc, use_unicode=True):
    docdicts = []
    for v in doc.itervalues():
        # check encoding of data unit
        if v.unit and use_unicode:
            data_unit = v.unit
        elif v.unit and not use_unicode:
            data_unit = v.unit.encode('utf8')
        else:
            data_unit = v.unit

        if isinstance(v.data, bool):
            vdata = str(v.data)
        elif isinstance(v.data, (datetime.date, datetime.time,
                                 datetime.datetime)):
            vdata = v.data.isoformat()
        else:
            vdata = v.data

        # append dictionary (one per value) to docdicts
        docdicts.append(
            {'Path': v.parent.parent.get_path(),
             'Section Name': v.parent.parent.name,
             'Section Definition': v.parent.parent.definition,
             'Section Type': v.parent.parent.type,
             'Property Name': v.parent.name,
             'Property Definition': v.parent.definition,
             'Value': vdata,
             'odML Data Type': v.dtype,
             'Data Unit':  data_unit,
             'Data Uncertainty': v.uncertainty})

    docdicts_sorted = sorted(docdicts, key=lambda k: k['Path'])

    return docdicts_sorted


def save_as_csv(docdicts, save_to):
    header = ['Path', 'Property Name', 'Value', 'Data Unit',
              'Data Uncertainty', 'odML Data Type', 'Property Definition']
    csvfile = open(save_to, 'wb')
    csvwriter = csv.DictWriter(csvfile, fieldnames=header, dialect='excel',
                               quoting=csv.QUOTE_NONNUMERIC)
    csvwriter.writeheader()
    for i, dnow in enumerate(docdicts):
        d = dict((k, '') for k in header)
        if i > 0:
            dold = docdicts[i-1]
            for k in header:
                if dnow[k] != dold[k]:
                    d[k] = dnow[k]
        else:
            d = dict((k, dnow[k]) for k in header)

        csvwriter.writerow(d)


def save_as_xls(docdicts, sheet_name, save_to):
    header = ['Path', 'Property Name', 'Value', 'Data Unit',
              'Data Uncertainty', 'odML Data Type', 'Property Definition']

    workbook = xlwt.Workbook()
    sheet = workbook.add_sheet(sheet_name)

    colwidth = {}
    for col, h in enumerate(header):
        colwidth[col] = [256 * (len(h)+1)]
        for d in docdicts:
            dtypes = (int, float, types.NoneType)
            content = str(d[h]) if isinstance(d[h], dtypes) else d[h]
            colwidth[col].append(256 * (len(content)+1))
        colwidth[col] = max(colwidth[col])

    style_1 = 'font: bold 1, color white; ' \
              'pattern: pattern solid, fore_colour Gray80;'
    style_2 = {1: 'font: bold 1, color white; '
                  'pattern: pattern solid, fore_colour gray80;',
               -1: 'font: bold 1; '
                   'pattern: pattern solid, fore_colour Gray25;'}
    style_3 = {1: 'pattern: pattern solid, fore_colour white;',
               -1: 'pattern: pattern solid, fore_colour Gray25;'}
    style_4 = 'pattern: pattern solid, fore_colour Red;'

    for col, h in enumerate(header):
        sheet.col(col).width = colwidth[col]
        style = xlwt.easyxf(style_1)
        sheet.write(0, col, h, style)
        for row, dnow in enumerate(docdicts):
            d = dict((k, v) for k, v in dnow.items())
            if row == 0:
                i = 1
                style_col0 = xlwt.easyxf(style_2[i*-1])
                style_rest = xlwt.easyxf(style_3[i*-1])
            else:
                dold = docdicts[row-1]
                if dnow['Property Name'] == dold['Property Name']:
                    d['Property Name'] = ''
                    d['Property Definition'] = ''
                if dnow['Path'] == dold['Path']:
                    d['Path'] = ''
                    style_col0 = style_col0
                    style_rest = style_rest
                else:
                    i = i*-1
                    style_col0 = xlwt.easyxf(style_2[i*-1])
                    style_rest = xlwt.easyxf(style_3[i*-1])

            if col in [0, 1]:
                sheet.write(row+1, col, d[h], style_col0)
            else:
                print(d[h])
                if h == 'Value' and d[h] in odml_default_value():
                    sheet.write(row+1, col, d[h], xlwt.easyxf(style_4))
                else:
                    sheet.write(row+1, col, d[h], style_rest)

    workbook.save(save_to)


def load_xls(load_from, sheet_name):
    # open excel sheet for reading
    output = None
    if os.path.exists(load_from):
        xlsbook = xlrd.open_workbook(load_from)
        sheet = xlsbook.sheet_by_name(sheet_name)

        header = sheet.row(0)

        val_attr = ['Value', 'Data Unit', 'Data Uncertainty', 'odML Data Type']

        output = []
        output_idx = 0
        for row_idx in xrange(1, sheet.nrows):
            row = sheet.row(row_idx)
            rowdict = dict((k.value, v.value) for (k, v) in zip(header, row))
            for k in val_attr:
                rowdict[k] = [rowdict[k]]

            if row_idx > 1:
                for k, v in rowdict.items():
                    if v == u'' and k not in val_attr:
                        rowdict[k] = output[output_idx-1][k]

                if rowdict['Property Name'] == \
                        output[output_idx-1]['Property Name']:
                    for k in val_attr:
                        output[output_idx-1][k].extend(rowdict[k])
                    continue
                else:
                    output.append(rowdict)
                    output_idx += 1

            else:
                output.append(rowdict)
                output_idx += 1

        for d in output:
            for i, dtype in enumerate(d['odML Data Type']):
                if d['Data Unit'][i] == '':
                    d['Data Unit'][i] = None
                if dtype == 'date':
                    f = '%Y-%m-%d'
                    d['Value'][i] = datetime.datetime.strptime(d['Value'][i],
                                                               f).date()
                elif dtype == 'time':
                    f = '%H:%M:%S'
                    d['Value'][i] = datetime.datetime.strptime(d['Value'][i],
                                                               f).time()
                elif dtype == 'datetime':
                    f = '%Y-%m-%d %H:%M:%S'
                    d['Value'][i] = datetime.datetime.strptime(d['Value'][i],
                                                               f)
                elif dtype == 'float':
                    d['Value'][i] = float(d['Value'][i])
                    d['Data Uncertainty'][i] = float(d['Data Uncertainty'][i])\
                        if d['Data Uncertainty'][i] != '' else None
                elif dtype == 'int':
                    d['Value'][i] = int(d['Value'][i])
                    d['Data Uncertainty'][i] = int(d['Data Uncertainty'][i])\
                        if d['Data Uncertainty'][i] != '' else None
                elif dtype == 'boolean':
                    if d['Value'][i] == 'True':
                        d['Value'][i] = True
                    else:
                        d['Value'][i] = False
                else:
                    d['Value'][i] = str(d['Value'][i])

    return output


def load_blackrock_electrodes_info(load_from):
    """
    Loads electrodes_info.txt file into a dictionary.

    Args:
        filename (str):
            Filename (absolute path) of electrodes_info.txt to load.

    Returns:
        (dict):
            key (int):
                ID of electrode used in nev and nsX files.
            values (dict):
                ConnectorAlignedID (key, str):
                    (value, int) connector aligned id of electrode
                BankID (key, str):
                    (value, int) connector bank id of electrode
                PinID (key, str):
                    (value, int) connector pin id of electrode
                GridID (key, str):
                    (value, int) array grid id of electrode
                Impedance (key, str):
                    (value, quantity) factory impedance of electrode
    """
    output = None
    if os.path.exists(load_from):
        csv_file = csv.DictReader(open(load_from, 'rb'), delimiter='\t')

        csv_dicts = []
        for line in csv_file:
            csv_dicts.append(dict((k, v) for k, v in line.iteritems()))

        output = {}
        for d in csv_dicts:
            elid = int(d['ElectrodeID'])
            if elid > 0:
                elid_str = 'Electrode_%03d' % elid
                output[elid_str] = {}
                for k in d.keys():
                    if k == 'ElectrodeID':
                        output[elid_str]['ID'] = {'data': int(d[k])}
                    elif k == 'ConnectorAlignedID':
                        output[elid_str]['IDca'] = {'data': int(d[k])}
                    elif k.startswith('Imp'):
                        k2, unit = k.split('_')
                        output[elid_str][k2] = {'data': float(d[k]),
                                                'unit': unit}
                    else:
                        output[elid_str][k] = {'data': int(d[k])}

    return output


def load_brain_aligned_elids(load_from):
    """
    Loads electrodes_info.txt file into a dictionary.

    Args:
        filename (str):
            Filename (absolute path) of electrodes_info.txt to load.

    Returns:
        (dict):
            key (int):
                ID of electrode used in nev and nsX files.
            values (dict):
                ConnectorAlignedID (key, str):
                    (value, int) connector aligned id of electrode
                BankID (key, str):
                    (value, int) connector bank id of electrode
                PinID (key, str):
                    (value, int) connector pin id of electrode
                GridID (key, str):
                    (value, int) array grid id of electrode
                Impedance (key, str):
                    (value, quantity) factory impedance of electrode
    """
    output = None
    if os.path.exists(load_from):
        csv_file = csv.DictReader(open(load_from, 'rb'), delimiter='\t')

        csv_dicts = []
        for line in csv_file:
            csv_dicts.append(dict((k, v) for k, v in line.iteritems()))

        output = {}
        for d in csv_dicts:
            elid = int(d['ElectrodeID'])
            if elid > 0:
                elid_str = 'Electrode_%03d' % elid
                output[elid_str] = {}
                output[elid_str]['ID'] = {'data': int(d['ElectrodeID'])}
                output[elid_str]['IDba'] = {'data': int(d['BrainAlignedID'])}

    return output


def load_autoimpedance(load_from, subsession_name):
    csv_file = csv.reader(open(load_from), quotechar='"')

#    output = list(np.ones_like(ConnectorAlignment) * -1)
#    for line in reader:
#        if line[0].startswith('*'):
#            continue
#        else:
#            if int(re.findall(r'\d+', line[0])[0]) < 97:
#                idx = ConnectorAlignment.index(int(re.findall(r'\d+',
#                                                   line[0])[0]))
#                output[idx] = int(re.findall(r'\d+', line[0])[1])

#    return output


def load_analog_events(load_from):
    # load results from force event detection mat file if they exist
    output = None
    if os.path.exists(load_from):
        evdet = scipy.io.loadmat(load_from, squeeze_me=True)

        # dict which help transform between mat file and python routines
        py2mat = {86: 1, 89: 2, 101: 3, 106: 4, 149: 5, 154: 6, 166: 7, 169: 8}
        mat2py = dict((v, k) for k, v in py2mat.iteritems())

        ev_idx = {0: 'OT', 1: 'HS', 2: 'OR', 3: 'BTB'}

        # get occurring real trial types (excluding tial type 0)
        trtys = [int(k[-1]) for k in evdet.keys() if k.startswith('trialIDs')]

        # define and fill analog event container
        output = {'trials': {}, 'sensors': {}}
        for tt in trtys:
            if hasattr(evdet['trialIDs_%i' % tt], '__iter__'):
                trids_tt = evdet['trialIDs_%i' % tt]
            else:
                trids_tt = [evdet['trialIDs_%i' % tt]]
            trids_tt = np.atleast_1d(trids_tt)
            chtridx_tt = evdet['checked_trials_%i' % tt]
            aev_tt = evdet['times_%i' % tt]
            for i, t in enumerate(trids_tt):
                output['trials'][t] = {'GripForceSignals': {},
                                       'DisplacementSignal': {}}
                output['trials'][t]['GripForceSignals'] = \
                    dict((ev_idx[k], v) for (k, v) in enumerate(aev_tt[i][:4]))
                output['trials'][t]['DisplacementSignal'] = \
                    dict((ev_idx[k], v) for (k, v) in enumerate(aev_tt[i][4:]))
                output['trials'][t]['trialtype'] = mat2py[tt]
                if i+1 in chtridx_tt:
                    output['trials'][t]['checked'] = True
                else:
                    output['trials'][t]['checked'] = False
            output['sensors'][mat2py[tt]] = \
                int(str(evdet['analyzed_sigs_%i' % tt])[-1])

    return output


def load_loadforce(load_from):
    # load results from force event detection mat file if they exist
    output = None
    if os.path.exists(load_from):
        lfdet = scipy.io.loadmat(load_from, squeeze_me=True)['alltr']
        output = {}
        for row in lfdet:
            if int(row[1]) == 0:
                output['Trial_%03i' % int(row[0])] = 'LF'
            elif int(row[1]) == 1:
                output['Trial_%03i' % int(row[0])] = 'HF'
            else:
                output['Trial_%03i' % int(row[0])] = 'n.d.'

    return output


def load_rejections(load_from):
    output = None
    filelist = glob.glob(load_from + '*.hdf5')
    if len(filelist) > 0:
        output = {}
        for filename in filelist:
            rej = h5py.load_h5(filename)
            output[rej['fsets']['filtername']] = \
                {'Highcut': rej['fsets']['highcut'],
                 'Lowcut': rej['fsets']['lowcut'],
                 'Order': rej['fsets']['order'],
                 'RejElectrodes': rej['elrej']['bad_electrodes']
                 if len(rej['elrej']['bad_electrodes']) > 0
                 else [-1],
                 'RejTrials': rej['trrej']['bad_trials']
                 if len(rej['trrej']['bad_trials']) > 0
                 else [-1]}

    return output


def load_unitids(load_from):
    output = None
    load_from_b = '.'.join([load_from.split('.')[0] + '-test', 'txt'])
    if os.path.exists(load_from):
        filename = load_from
    elif os.path.exists(load_from_b):
        filename = load_from_b
    else:
        filename = None

    if filename:
        csv_file = csv.reader(open(filename), delimiter='\t')

        output = {}
        elids = []
        for row in csv_file:
            elid_str = row[0]
            if elid_str != 'NaN':
                if int(elid_str) in elids:
                    raise ValueError("Unitid txt-file is corrupt! Electrode "
                                     "ID %s is entered twice!" % elid_str)
                elids.append(int(elid_str))
                elid = 'Electrode_%03d' % int(elid_str)
                if int(row[2]) > 0:
                    if int(row[1]) + 1 != int(row[2]):
                        raise ValueError("Unitid txt-file is corrupt! %s: MUA "
                                         "ID must be one larger than the "
                                         "largest SUA ID" % elid)
                output[elid] = {'SUAIDs': range(1, int(row[1])+1)
                                if int(row[1]) > 0
                                else [-1],
                                'MUAID': int(row[2])
                                if int(row[2]) > 0
                                else -1}
        filename = os.path.basename(filename)
    else:
        filename = None

    return output, filename


def save_templates(save_to):
    filename = 'blackrock_template'
    doc = blackrock.blackrock_tree(connector_style='CerePort', grid_no=1)
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata blackrock', save_to + filename + '.xls')

    filename = 'blackrock_info_template'
    doc = blackrock.blackrock_tree(connector_style='CerePort', grid_no=1)
    doc.remove(doc['Headstage'])
    parent = doc['Cerebus']
    parent.remove(parent['NeuralSignalAmplifier'])
    parent = doc['Cerebus']['NeuralSignalProcessor']
    parent.remove(parent['AnalogSignals'])
    parent.remove(parent['DigitalSignals'])
    parent.remove(parent['NeuralSignals'])
    parent = doc['Cerebus']['NeuralSignalProcessor']['SpikeDetection']
    parent.remove(parent.properties['ThresholdType'])
    parent.remove(parent.properties['SavedIn'])
    parent.remove(parent['Events'])
    parent.remove(parent['Filter'])
    parent.remove(parent['Rejections'])
    parent = doc['UtahArray']['Array']
    parent.remove(parent['Electrode_XXX'])
    parent.remove(parent.properties['UsedElectrodeCount'])
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata blackrock', save_to + filename + '.xls')

    filename = 'project_template'
    doc = project.rg_project_tree()
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata project', save_to + filename + '.xls')

    filename = 'subject_template'
    doc = subject.rg_subject_tree()
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata subject', save_to + filename + '.xls')

    filename = 'setup_template'
    doc = setup.rg_setup_tree()
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata setup', save_to + filename + '.xls')

#    filename = 'preprocessing_template'
#    doc = preprocessing.rg_preprocessing_tree()
#    docdicts = write_odml2dict(doc, use_unicode=False)
#    save_as_csv(docdicts, save_to + filename + '.csv')
#    docdicts = write_odml2dict(doc, use_unicode=True)
#    save_as_xls(docdicts, 'metadata', save_to + filename + '.xls')

    filename = 'subsession_template'
    doc = subsession.rg_subsession_tree()
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata subsession', save_to + filename + '.xls')

    filename = 'subsession_template_OneCue'
    doc = subsession.rg_subsession_tree(tasktype='OneCue')
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata subsession', save_to + filename + '.xls')

    filename = 'subsession_template_TwoCues'
    doc = subsession.rg_subsession_tree(tasktype='TwoCues')
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata subsession', save_to + filename + '.xls')

    filename = 'subsession_template_Mapping'
    doc = subsession.rg_subsession_tree(tasktype='Mapping')
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata subsession', save_to + filename + '.xls')

    filename = 'subsession_template_Observation'
    doc = subsession.rg_subsession_tree(tasktype='Observation')
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata subsession', save_to + filename + '.xls')

    filename = 'subsession_template_Sleep'
    doc = subsession.rg_subsession_tree(tasktype='Sleep')
    rgodml.odml_io.save_odml(doc, save_to + filename + '.odml')
    docdicts = write_odml2dict(doc, use_unicode=False)
    save_as_csv(docdicts, save_to + filename + '.csv')
    docdicts = write_odml2dict(doc, use_unicode=True)
    save_as_xls(docdicts, 'metadata subsession', save_to + filename + '.xls')