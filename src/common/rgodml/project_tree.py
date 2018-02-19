# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:28:19 2014

@author: zehl
"""
import odml
from rgodml.odml_io import default_value


def rg_collaborators():
    """
    Returns a list of collaborators of the reach-to-grasp project

    Returns:
        (list of strings)
            Affiliations of involved collaborators
    """

    col = ['Inst. of Neuroscience and Medicine (INM-6) and Inst. for '
           'Advanced Simulation (IAS-6), Juelich Research Centre & JARA, '
           'Juelich, Germany',
           'Inst. de Neurosciences de la Timone (INT), UMR 7289, CNRS - '
           'Aix Marseille Univ., Marseille, France']

    return col


def rg_project_tree(tasktype='TwoCues'):
    """
    Builds an odML Document for reach-to-grasp project metadata

    Args:
        used_taskdesign (string, optional):
            States the task design used in this subsession

    Returns:
        (odml.doc.BaseDocument)
            odML Document with default reach-to-grasp project metadata

    Note:
        Valid used_taskdesign: OneCue, TwoCues, Observation, Sleep, Mapping
    """

    # CREATE A DOCUMENT
    doc = odml.Document(version='0.1')

    # CREATE A MAIN SECTION
    doc.append(odml.Section(name='Project',
                            type='project',
                            definition='Information on the experimental '
                                       'project'))

    # SET PARENT NODE
    parent = doc['Project']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='TaskDesigns',
                               type='project',
                               definition='Information on tasks'))

    parent.append(odml.Section(name='TrialTypeCodes',
                               type='codes',
                               definition='Information on trial type codes'))

    parent.append(odml.Section(name='ConditionCodes',
                               type='codes',
                               definition='Information on condition codes '
                                          '(trial type combination codes)'))

    parent.append(odml.Section(name='PerformanceCodes',
                               type='codes',
                               definition='Information on trial performance '
                                          'codes'))

    parent.append(odml.Section(name='TrialSettings',
                               type='setup/control',
                               definition='Information on the programmed '
                                          'trial settings'))

    parent.append(odml.Section(name='TargetObjectSettings',
                               type='setup/control',
                               definition='Information on the settings of the '
                                          'target object'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Name',
                                value=odml.Value(data='reach-to-grasp',
                                                 dtype=odml.DType.string),
                                definition='Name of project'))

    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='electrophysiology',
                                                 dtype=odml.DType.string),
                                definition='Type of project'))

    parent.append(odml.Property(name='Subtype',
                                value=odml.Value(data='motor behavior',
                                                 dtype=odml.DType.string),
                                definition='Name of the project'))

    parent.append(odml.Property(name='Description',
                                value=default_value(odml.DType.url),
                                definition='Project description URL'))

    parent.append(odml.Property(name='Collaborators',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in rg_collaborators()],
                                definition='Affilitations of involved '
                                           'collaborators'))

    # SET PARENT NODE
    parent = doc['Project']['TaskDesigns']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Designs',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['OneCue', 'TwoCues',
                                                    'Observation', 'Sleep',
                                                    'Mapping']],
                                definition='Names of possible task designs'))

    parent.append(odml.Property(name='UsedDesign',
                                value=odml.Value(data=tasktype,
                                                 dtype=odml.DType.string)
                                if tasktype
                                else default_value(odml.DType.string),
                                definition='Name of used task design'))

    parent.append(odml.Property(name='Description',
                                value=default_value(odml.DType.url),
                                definition='Task designs descriptions URL'))

    # SET PARENT NODE
    parent = doc['Project']['TrialTypeCodes']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='SGLF',
                                value=odml.Value(data=89,
                                                 dtype=odml.DType.int),
                                definition='Side grip / low force code'))

    parent.append(odml.Property(name='SGHF',
                                value=odml.Value(data=86,
                                                 dtype=odml.DType.int),
                                definition='Side grip / high force code'))

    parent.append(odml.Property(name='PGLF',
                                value=odml.Value(data=169,
                                                 dtype=odml.DType.int),
                                definition='Precision grip / low force code'))

    parent.append(odml.Property(name='PGHF',
                                value=odml.Value(data=166,
                                                 dtype=odml.DType.int),
                                definition='Precision grip / high force code'))

    parent.append(odml.Property(name='LFSG',
                                value=odml.Value(data=149,
                                                 dtype=odml.DType.int),
                                definition='Low force / side grip code'))

    parent.append(odml.Property(name='LFPG',
                                value=odml.Value(data=154,
                                                 dtype=odml.DType.int),
                                definition='Low force / precision grip code'))

    parent.append(odml.Property(name='HFSG',
                                value=odml.Value(data=101,
                                                 dtype=odml.DType.int),
                                definition='High force / side grip code'))

    parent.append(odml.Property(name='HFPG',
                                value=odml.Value(data=106,
                                                 dtype=odml.DType.int),
                                definition='High force / precision grip code'))

    parent.append(odml.Property(name='SGSG',
                                value=odml.Value(data=85,
                                                 dtype=odml.DType.int),
                                definition='Side grip / side grip code'))

    parent.append(odml.Property(name='PGPG',
                                value=odml.Value(data=170,
                                                 dtype=odml.DType.int),
                                definition='Precision grip / precision grip '
                                           'code'))

    parent.append(odml.Property(name='n.d.',
                                value=odml.Value(data=0,
                                                 dtype=odml.DType.int),
                                definition='Code for not detectable trial '
                                           'type'))

    # SET PARENT NODE
    parent = doc['Project']['ConditionCodes']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='cnd_1',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGLF', 'SGHF',
                                                    'PGLF', 'PGHF']],
                                definition='Grip first trials, all grip and '
                                           'force types'))

    parent.append(odml.Property(name='cnd_2',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGLF', 'SGHF',
                                                    'PGLF', 'PGHF']],
                                definition='Force first trials, all force and '
                                           'grip types'))

    parent.append(odml.Property(name='cnd_11',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGLF', 'PGLF']],
                                definition='Grip first trials, all grip '
                                           'types, but only low force'))

    parent.append(odml.Property(name='cnd_12',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGHF', 'PGHF']],
                                definition='Grip first trials, all grip '
                                           'types, but only high force'))

    parent.append(odml.Property(name='cnd_13',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGHF', 'SGLF']],
                                definition='Grip first trials, only side '
                                           'grip, but all force types'))

    parent.append(odml.Property(name='cnd_14',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['PGHF', 'PGLF']],
                                definition='Grip first trials, only precision '
                                           'grip, but all force types'))

    parent.append(odml.Property(name='cnd_21',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['LFSG', 'LFPG']],
                                definition='Force first trials, only low '
                                           'force, but all grip types'))

    parent.append(odml.Property(name='cnd_22',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['HFSG', 'HFPG']],
                                definition='Force first trials, only high '
                                           'force, but all grip types'))

    parent.append(odml.Property(name='cnd_23',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['HFSG', 'LFSG']],
                                definition='Force first trials, all force '
                                           'types, but only side grip'))

    parent.append(odml.Property(name='cnd_24',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['HFPG', 'LFPG']],
                                definition='Force first trials, all force '
                                           'types, but only side grip'))

    parent.append(odml.Property(name='cnd_131',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGLF']],
                                definition='Grip first trials, only side '
                                           'grip, only low force'))

    parent.append(odml.Property(name='cnd_132',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGHF']],
                                definition='Grip first trials, only side '
                                           'grip, only high force'))

    parent.append(odml.Property(name='cnd_133',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['SGSG']],
                                definition='Grip first trials, only side '
                                           'grip, force needs to be guessed'))

    parent.append(odml.Property(name='cnd_141',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['PGLF']],
                                definition='Grip first trials, only precision '
                                           'grip, only low force'))

    parent.append(odml.Property(name='cnd_142',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['PGHF']],
                                definition='Grip first trials, only precision '
                                           'grip, only high force'))

    parent.append(odml.Property(name='cnd_144',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['PGPG']],
                                definition='Grip first trials, only precision '
                                           'grip, force needs to be guessed'))

    parent.append(odml.Property(name='cnd_213',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['LFSG']],
                                definition='Force first trials, only low '
                                           'force, only side grip'))

    parent.append(odml.Property(name='cnd_214',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['LFPG']],
                                definition='Force first trials, only low '
                                           'force, only precision grip'))

    parent.append(odml.Property(name='cnd_223',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['HFSG']],
                                definition='Force first trials, only high '
                                           'force, only side grip'))

    parent.append(odml.Property(name='cnd_224',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['HFPG']],
                                definition='Force first trials, only high '
                                           'force, only precision grip'))

    # SET PARENT NODE
    parent = doc['Project']['PerformanceCodes']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='pc_127',
                                value=odml.Value(data='correct trial',
                                                 dtype=odml.DType.string),
                                definition='Performance code of a correctly '
                                           'completed trial'))

    parent.append(odml.Property(name='pc_64',
                                value=odml.Value(data='SR < FP-ON',
                                                 dtype=odml.DType.string),
                                definition='Performance code of a trial where '
                                           'switch is released too early '
                                           '(before fixation point)'))

    parent.append(odml.Property(name='pc_96',
                                value=odml.Value(data='SR < CUE-ON',
                                                 dtype=odml.DType.string),
                                definition='Performance code of a trial where '
                                           'switch is released too early '
                                           '(before delay cue is turned on)'))

    parent.append(odml.Property(name='pc_112',
                                value=odml.Value(data='SR < CUE-OFF',
                                                 dtype=odml.DType.string),
                                definition='Performance code of a trial where '
                                           'switch is released too early '
                                           '(before delay cue is turned off)'))

    parent.append(odml.Property(name='pc_120',
                                value=odml.Value(data='SR < GO-ON',
                                                 dtype=odml.DType.string),
                                definition='Performance code of a trial where '
                                           'switch is released too early '
                                           '(before GO cue turns on)'))

    parent.append(odml.Property(name='pc_124',
                                value=odml.Value(data='no SR',
                                                 dtype=odml.DType.string),
                                definition='Performance code of a trial where '
                                           'switch is not released'))

    parent.append(odml.Property(name='pc_126',
                                value=odml.Value(data='grip error',
                                                 dtype=odml.DType.string),
                                definition='Performance code of a trial where '
                                           'the wrong grip type is used'))

    # SET PARENT NODE
    parent = doc['Project']['TrialSettings']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='DigitalEvents',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['TS', 'FP-ON', 'CUE-ON',
                                                    'CUE-OFF', 'GO-ON', 'SR',
                                                    'RW', 'GO-OFF']],
                                definition='Digital recorded trial events'))

    parent.append(odml.Property(name='AnalogEvents',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['OT', 'HS', 'OR', 'BTB']],
                                definition='Analog recorded trial events'))

    parent.append(odml.Property(name='TSP',
                                value=odml.Value(data=400.0,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Duration of trial start period '
                                           '(period between TS and FP)'))

    parent.append(odml.Property(name='FPP',
                                value=odml.Value(data=400.0,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Duration of fixation point period '
                                           '(period between FP and CUE-ON)'))

    parent.append(odml.Property(name='CUE',
                                value=odml.Value(data=300.0,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Duration of CUE period (period '
                                           'between CUE-ON and CUE-OFF'))

    parent.append(odml.Property(name='DELAY',
                                value=odml.Value(data=1000.0,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Duration of delay period (period '
                                           'between CUE-OFF and GO)'))

    parent.append(odml.Property(name='HT',
                                value=odml.Value(data=500.0,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Duration of holding time (period '
                                           'between HS and RW)'))

    parent.append(odml.Property(name='MTLimit',
                                value=odml.Value(data=700.0,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Time limit for movement (period '
                                           'between GO-ON and HS)'))

    parent.append(odml.Property(name='RTLimit',
                                value=odml.Value(data=1000.0,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Time limit for reaction time '
                                           '(period between GO-ON and SR)'))

    parent.append(odml.Property(name='PTLimits',
                                value=default_value(odml.DType.float,
                                                    unit='ms'),
                                definition='Min and max time to pull the '
                                           'object into the holding position '
                                           '(period between OT and HS)'))

    # SET PARENT NODE
    parent = doc['Project']['TargetObjectSettings']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='HoldPosRange',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.int,
                                                  unit='mm')
                                       for data in [4.0, 14.0]],
                                definition='Distance window for holding the '
                                           'object in place'))

    parent.append(odml.Property(name='LowForceWeight',
                                value=odml.Value(data=100.0,
                                                 dtype=odml.DType.float,
                                                 unit='g'),
                                definition='Weight of object in low force '
                                           'modus'))

    parent.append(odml.Property(name='HighForceWeight',
                                value=odml.Value(data=200.0,
                                                 dtype=odml.DType.float,
                                                 unit='g'),
                                definition='Weight of object in high force '
                                           'modus'))

    return doc
