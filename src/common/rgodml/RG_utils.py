# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:32:12 2014

@author: zehl
"""
import os
import re
import csv
import glob

import numpy as np

from util import hdf5_wrapper
from util import excel_reader


def get_autoimpedance(directory, subsession, ConnectorAlignment):
    avail_autoimps = np.sort(get_filenames(directory, "*impedance.txt"))

    subses_date = subsession[:7]

    filename = None
    for autoimp_date in [i[:7] for i in avail_autoimps]:
        if autoimp_date <= subses_date:
            filename = os.path.join(directory, autoimp_date + "_impedance.txt")

    if filename != None:
        output = {'Filename': filename,
                  'Impedances': list(read_blr_autoimpedance(filename, ConnectorAlignment))}
    else:
        output = {'Filename': '-',
                  'Impedances': 0.0}

    return output

def read_blr_electrode_info(filename):
    reader = csv.DictReader(open(filename), delimiter='\t', quotechar='"')

    output = {}
    for linedict in reader:
        for key, value in linedict.iteritems():
            output.setdefault(key, []).append(value)

    for key in output.keys():
        try:
            if key in ['Impedances_MOhm']:
                output[key] = [float(i) for i in output[key]]
            else:
                output[key] = [int(i) for i in output[key]]
        except:
            pass

    return output

def read_brain_alignment(filename):
    reader = csv.DictReader(open(filename), delimiter='\t', quotechar='"')

    output = {}
    for linedict in reader:
        for key, value in linedict.iteritems():
            output.setdefault(key, []).append(int(value))

    for key in output.keys():
        try:
            output[key] = [float(i) for i in output[key]]
        except:
            pass

    return output

def read_blr_autoimpedance(filename, ConnectorAlignment):
    reader = csv.reader(open(filename), quotechar='"')

    output = list(np.ones_like(ConnectorAlignment) * -1)
    for line in reader:
        if line[0].startswith('*'):
            continue
        else:
            if int(re.findall(r'\d+', line[0])[0]) < 97:
                idx = ConnectorAlignment.index(int(re.findall(r'\d+', line[0])[0]))
                output[idx] = int(re.findall(r'\d+', line[0])[1])

    return output

def read_small_xls(filename):
    output = {}

    md = excel_reader.xls_sheet_reader(filename, 'Sheet1', 2, '', 1)[0]
    for dk in np.sort(md.keys()):
        if str(md[dk]['Value']).strip() != '':
            dtype = md[dk]['dtype']
            if dtype == 'date':
                value = excel_reader.datefloat2date(md[dk]['Value'])
            elif dtype in ['int', 'float']:
                if isinstance(md[dk]['Value'], (int, float)):
                    value = md[dk]['Value']
                else:
                    value = excel_reader.string2list(md[dk]['Value'], dtype)
            else:
                value = md[dk]['Value']
            output[dk] = value

    return output

def read_subsession_info(filename, subsession):
    output = {}

    md = excel_reader.xls_sheet_reader(filename, 'Sheet1', 2, '-', 7)[0][subsession]

    md = hdf5_wrapper.flatten(md, separator='/')

    output = {}
    for k in md.keys():
        if k == 'Behavior/Comment':
            output['/Experiment/' + k] = md[k]
        elif type(md[k]) == float:
            output['/Experiment/' + k] = md[k]
        elif str(md[k]).strip() == '':
            output['/Experiment/' + k] = 0.0
        elif k in ['Behavior/Noisy', 'Behavior/Task/OneCue', 'Behavior/Task/TwoCues',
                   'Behavior/Task/Observation', 'Behavior/Settings/CueTask/Standard',
                   'Behavior/Settings/CueTask/TrialTypeOrder/Block/Grip',
                   'Behavior/Settings/CueTask/TrialTypeOrder/Block/Force',
                   'Behavior/Settings/CueTask/TrialTypeOrder/Random/Grip',
                   'Behavior/Settings/CueTask/TrialTypeOrder/Random/Force',
                   'Behavior/Task/Sleep', 'Behavior/Task/Mapping']:
#            if md[k].strip() == 'T':
            output['/Experiment/' + k] = 'True' if md[k].strip() == 'T' else 'False'
#            elif md[k].strip() == 'F':
#                output['/Experiment/'+k] = 'False'
        else:
            output['/Experiment/' + k] = [float(i) for i in re.findall(r'\d+.\d+|\d+', md[k])]

    return output


def get_filenames(directory, filename):

    filepaths = glob.glob(os.path.join(directory, filename))
    output = [os.path.basename(i) for i in filepaths]

    return output

def get_spike_amplitude_thresholds(subsession_object):
    el_pa = subsession_object.parameters_electrodes

    idx = [idx for idx, pe in enumerate(el_pa) if pe.has_key('ConnectorID') and pe['ConnectorID'] <= 3]
    blr_id = [el_pa[i]['ConnectorID'] * el_pa[i]['ConnectorPin'] for i in idx]

    output = list(np.zeros_like(idx))
    for i, e in zip(idx, blr_id):
        output[e - 1] = el_pa[i]['AmpThresholdLo']

    return output
