# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 12:23:54 2015

@author: zehl
"""
import time
t1 = time.time()

import odml
import datetime

import rg.rgio

# import template scripts for odML branches
import rgodml.setup_tree
import rgodml.subject_tree
import rgodml.project_tree
import rgodml.blackrock_tree
import rgodml.subsession_tree
import rgodml.preprocessing_tree

# import rg odml utils
import rgodml.metadata_io
import rgodml.odml_utils
import rgodml.odml_io

# =============================================================================
# DEFINE MONKEY AND PATHS
# =============================================================================
monkey = 'Nikos2'
subsession_name = 'i140725-002'
#data_dir = '/media/zehl/DATA/datasets/DataGrasp/Data%s/' % monkey
#mdata_dir = '/media/zehl/DATA/datasets/DataGrasp/MetaData%s/' % monkey
data_dir = '/riou/work/comco_congloue/DataGrasp/Data%s/' % monkey
mdata_dir = '/riou/work/comco_congloue/DataGrasp/MetaData%s/' % monkey
odml_dir = mdata_dir + '/odMLfiles/'

# =============================================================================
# LOAD METADATA SOURCES
# =============================================================================
try:
    subsobj = rg.rgio.ReachGraspIO(data_dir + subsession_name)
except:
    subsobj = None

load_from = mdata_dir + '/source/subsessions/' + subsession_name + '.xls'
subsxls = rgodml.metadata_io.load_xls(load_from, 'metadata subsession')

load_from = mdata_dir + '/source/monkey/monkey_info.xls'
monkeyxls = rgodml.metadata_io.load_xls(load_from, 'metadata subject')

load_from = mdata_dir + '/source/blackrock/blackrock_info.xls'
brxls = rgodml.metadata_io.load_xls(load_from, 'metadata blackrock')

load_from = mdata_dir + '/source/blackrock/electrodes_info.txt'
elinfo = rgodml.metadata_io.load_blackrock_electrodes_info(load_from)

load_from = mdata_dir + '/source/monkey/brain_aligned_elids.txt'
ba_elids = rgodml.metadata_io.load_brain_aligned_elids(load_from)

load_from = mdata_dir + '/source/event_detect/fe_' + subsession_name + '_m.mat'
aev = rgodml.metadata_io.load_analog_events(load_from)

load_from = mdata_dir + '/source/load_detect/lf_' + subsession_name + '.mat'
alf = rgodml.metadata_io.load_loadforce(load_from)

load_from = mdata_dir + '/source/rejections/data/' + subsession_name
lfprej = rgodml.metadata_io.load_rejections(load_from)

load_from = data_dir + subsession_name + '.txt'
unitids, unitids_filename = rgodml.metadata_io.load_unitids(load_from)

# =============================================================================
# LOAD (SPECIFIED) TEMPLATES
# =============================================================================

# SUBJECT TEMPLATE
subject = rgodml.subject_tree.rg_subject_tree()

# SETUP TEMPLATE
# location = 'old' for Lilou and Tanya
# location = 'new' for Tanya2, Nikos, and Nikos2
setup = rgodml.setup_tree.rg_setup_tree(location='old')

# PROJECT TEMPLATE
# tasktype in subsxls specifies the template
if subsxls:
    tasktype = [d['Value'][0] for d in subsxls if
                d['Path'] == '/Subsession' and
                d['Property Name'] == 'TaskType'][0]
else:
    tasktype = None
project = rgodml.project_tree.rg_project_tree(tasktype=tasktype)

# BLACKROCK TEMPLATE
# electrode_ids, avail_nsx and analogI_chIDs in subsobj specify the template
if subsobj:
    electrode_ids = [e for e in subsobj.channel_id_nev if e <= 96]
    avail_nsx = subsobj.nsx_avail
    if 2 in avail_nsx:
        analogI_chIDs = [i for i in subsobj.channel_id_nsx[2] if i > 128]
    else:
        analogI_chIDs = []
else:
    electrode_ids = []
    avail_nsx = []
    analogI_chIDs = []
if unitids:
    spike_sorted = True
else:
    spike_sorted = False
# set parameters to specify the template
# grid_no = 1 for all subjects
# connector_style = 'CerePort' for Lilou, Tanya, and Tanya2
# connector_style = ??? for Nikos and Nikos2
if subsxls:
    used_headstage = [d['Value'][0] for d in subsxls if
                      d['Path'] == '/Headstage' and
                      d['Property Name'] == 'Model'][0]
else:
    used_headstage = None
blackrock = \
    rgodml.blackrock_tree.blackrock_tree(used_DIOports=['ExpDigitalInput'],
                                         connector_style='CerePort',
                                         grid_no=1,
                                         electrode_ids=electrode_ids,
                                         avail_nsx=avail_nsx,
                                         analogI_chIDs=analogI_chIDs,
                                         used_headstage=used_headstage,
                                         spike_sorted=spike_sorted)

# SUBSESSION TEMPLATE
# tasktype in subsxls and trial ids in subsobj specify template settings
if subsobj:
    trialids = subsobj.get_trial_ids()
else:
    trialids = []
subsession = rgodml.subsession_tree.rg_subsession_tree(tasktype=tasktype,
                                                       trialids=trialids)

# LFP REJECTION TEMPLATE
if lfprej:
    rej_names = lfprej.keys()
else:
    rej_names = []
rejections = rgodml.preprocessing_tree.rej_tree(rej_names=rej_names)
# =============================================================================
# CREATE UNIFIED ODML DOCUMENT FOR SUBSESSION
# =============================================================================

# Define odml document properties
doc = odml.Document(author='Lyuba Zehl',
                    version='0.1',
                    date=datetime.date.today())

# Append all branches to odml document
for s in setup.sections:
    doc.append(s.clone())
for s in subject.sections:
    doc.append(s.clone())
for s in project.sections:
    doc.append(s.clone())
for s in blackrock.sections:
    doc.append(s.clone())
for s in subsession.sections:
    doc.append(s.clone())
for s in rejections.sections:
    doc.append(s.clone())

# =============================================================================
# FILL UNIFIED ODML DOCUMENT
# =============================================================================

# with metadata from subsession_name.xls
if subsxls:
    source_filename = subsession_name + '.xls'
    rgodml.odml_utils.fill_odml_with_xls(doc, subsxls, source_filename)

# with metadata from monkey_info.xls
if monkeyxls:
    source_filename = 'monkey_info.xls'
    rgodml.odml_utils.fill_odml_with_xls(doc, monkeyxls, source_filename)

# with metadata from blackrock_info.xls
if brxls:
    source_filename = 'blackrock_info.xls'
    rgodml.odml_utils.fill_odml_with_xls(doc, brxls, source_filename)

# with metadata from electrodes_info.txt
if elinfo and len(electrode_ids) > 0:
    source_filename = 'electrodes_info.txt'
    rgodml.odml_utils.fill_odml_with_elinfo(doc, elinfo, source_filename)

# with metadata from brain_aligned_elids.txt
if ba_elids and len(electrode_ids) > 0:
    source_filename = 'brain_aligned_elids.txt'
    rgodml.odml_utils.fill_odml_with_ba_elids(doc, ba_elids, source_filename)

# with metadata from analog event detection mat file
if aev:
    source_filename = 'fe_' + subsession_name + '_m.mat'
    rgodml.odml_utils.fill_odml_with_aev(doc, aev, source_filename)

# with metadata from analog force load detection mat file
if alf:
    source_filename = 'lf_' + subsession_name + '.mat'
    rgodml.odml_utils.fill_odml_with_alf(doc, alf, source_filename)

# with metadata from subsession_name.nev
if subsobj:
    source_filename = subsession_name + '.nev'
    rgodml.odml_utils.fill_odml_with_subsobj(doc, subsobj, source_filename)

# with metadata from spikesorting (unitids)
if unitids:
    rgodml.odml_utils.fill_odml_with_unitids(doc, unitids, unitids_filename)

# with metadata for periods
if subsobj and aev:
    rgodml.odml_utils.fill_odml_with_periods(doc, subsobj, aev)

# with metadata for lfp rejections
if lfprej:
    rgodml.odml_utils.fill_odml_with_lfprej(doc, lfprej, subsession_name)

# =============================================================================
# Save odml Document
# =============================================================================

save_to = odml_dir + subsession_name + '.odml'
rgodml.odml_io.save_odml(doc, save_to)

t2 = time.time()
print (t2-t1), 's'
