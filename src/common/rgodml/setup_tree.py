# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:38:50 2014

@author: zehl
"""
import odml
from rgodml.odml_io import default_value


def get_location(location):
    if location == 'old':
        loc = odml.Value(data='Inst. de Neurosciences Cognitives de la '
                              'Mediterranee (INCM), GLM, CNRS - Aix Marseille '
                              'Univ., Marseille, France',
                         dtype=odml.DType.string)
    elif location == 'new':
        loc = odml.Value(data='Inst. de Neurosciences de la Timone (INT), '
                              'UMR 7289, CNRS - Aix Marseille Univ., '
                              'Marseille, France',
                         dtype=odml.DType.string)
    else:
        loc = default_value(odml.DType.url)

    return loc


def rg_setup_tree(location='-'):
    """
    Builds an odML Document for reach-to-grasp setup metadata

    Returns:
        (odml.doc.BaseDocument)
            odML Document with default reach-to-grasp setup metadata
    """
    # CREATE A DOCUMENT
    doc = odml.Document(version='0.1')

    # APPEND MAIN SECTION
    doc.append(odml.Section(name='Setup',
                            type='setup',
                            definition='Information on the experimental '
                                       'setup'))

    # SET PARENT NODE
    parent = doc['Setup']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='Apparatus',
                               type='setup',
                               definition='Information on the experimental '
                                          'apparatus'))

    parent.append(odml.Section(name='ControlComputer',
                               type='setup/software',
                               definition='Information on the experimental '
                                          'control'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Location',
                                value=get_location(location),
                                definition='Location of the setup'))

    parent.append(odml.Property(name='Description',
                                value=default_value(odml.DType.url),
                                definition='Setup description URL'))

    parent.append(odml.Property(name='DAQSystem',
                                value=odml.Value(data='Cerebus',
                                                 dtype=odml.DType.string),
                                definition='DAQ system used'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='StartSwitch',
                               type='setup',
                               definition='Information on the start switch'))

    parent.append(odml.Section(name='TargetObject',
                               type='setup',
                               definition='Information on the target object'))

    parent.append(odml.Section(name='CueSystem',
                               type='setup',
                               definition='Information on the cue system'))

    parent.append(odml.Section(name='RewardSystem',
                               type='setup',
                               definition='Information on the reward system'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Manufacturer',
                                value=odml.Value(data='NI A/D Card (ref?)',
                                                 dtype=odml.DType.string),
                                definition='Owner of the setup'))

    parent.append(odml.Property(name='Creator',
                                value=odml.Value(data='Thomas Brochier',
                                                 dtype=odml.DType.string),
                                definition='Full name of person who created '
                                           'the setup'))

    parent.append(odml.Property(name='Maintainer',
                                value=odml.Value(data='Thomas Brochier',
                                                 dtype=odml.DType.string),
                                definition='Full name of person who maintains '
                                           'the setup'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['StartSwitch']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Function',
                                value=odml.Value(data='activates start trial',
                                                 dtype=odml.DType.string),
                                definition='Function of apparatus'))

    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='table switch',
                                                 dtype=odml.DType.string),
                                definition='Type of apparatus'))

    parent.append(odml.Property(name='SignalType',
                                value=odml.Value(data='analog',
                                                 dtype=odml.DType.string),
                                definition='Signal type of apparatus (digital '
                                           'or analog)'))

    parent.append(odml.Property(name='Events',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['TS', 'SR']],
                                definition='Names of events extracted from '
                                           'this signal'))

    parent.append(odml.Property(name='SamplingRate',
                                value=odml.Value(data=30000,
                                                 dtype=odml.DType.int,
                                                 unit='samples/sec'),
                                definition='Sampling rate of data'))

    parent.append(odml.Property(name='ConnectedTo',
                                value=odml.Value(data='NeuralSignalProcessor',
                                                 dtype=odml.DType.string),
                                definition='Target device of apparatus'))

    parent.append(odml.Property(name='MonitoredBy',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Control software of apparatus'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='nev',
                                                 dtype=odml.DType.string),
                                definition='File in which apparatus signals '
                                           'are saved'))

    parent.append(odml.Property(name='Length',
                                value=default_value(odml.DType.float,
                                                    unit='cm'),
                                definition='Length of apparatus'))

    parent.append(odml.Property(name='Width',
                                value=default_value(odml.DType.float,
                                                    unit='cm'),
                                definition='Length of apparatus'))

    parent.append(odml.Property(name='Height',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Length of apparatus'))

    parent.append(odml.Property(name='PosX',
                                value=odml.Value(data=5.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='X position of apparatus'))

    parent.append(odml.Property(name='PosY',
                                value=odml.Value(data=0.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='Y position of apparatus'))

    parent.append(odml.Property(name='PosZ',
                                value=odml.Value(data=0.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='Z position of apparatus'))

    parent.append(odml.Property(name='PosOrigin',
                                value=odml.Value(data='midline, waist-level',
                                                 dtype=odml.DType.string),
                                definition='Origin for xyz position of '
                                           'apparatus'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['TargetObject']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='FSRSensor',
                               type='setup/daq',
                               definition='Information on the FSR sensors'))

    parent.append(odml.Section(name='HESensor',
                               type='setup/daq',
                               definition='Information on the hall-effect '
                                          'sensor'))

    parent.append(odml.Section(name='LoadForce',
                               type='setup/daq',
                               definition='Information on the load force '
                                          'signal'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Function',
                                value=odml.Value(data='object to grasp, pull '
                                                      'and hold',
                                                 dtype=odml.DType.string),
                                definition='Function of apparatus'))

    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='stainless steel '
                                                      'parallelpiped',
                                                 dtype=odml.DType.string),
                                definition='Type of apparatus'))

    parent.append(odml.Property(name='Length',
                                value=odml.Value(data=40.0,
                                                 dtype=odml.DType.float,
                                                 unit='mm'),
                                definition='Length of apparatus'))

    parent.append(odml.Property(name='Width',
                                value=odml.Value(data=16.0,
                                                 dtype=odml.DType.float,
                                                 unit='mm'),
                                definition='Width of apparatus'))

    parent.append(odml.Property(name='Height',
                                value=odml.Value(data=10.0,
                                                 dtype=odml.DType.float,
                                                 unit='mm'),
                                definition='Height of apparatus'))

    parent.append(odml.Property(name='AttachedTo',
                                value=odml.Value(data='anterior end of a low- '
                                                      'friction horizontal '
                                                      'shuttle',
                                                 dtype=odml.DType.string),
                                definition='Device apparatus is attached to'))

    parent.append(odml.Property(name='Rotation',
                                value=odml.Value(data=45.0,
                                                 dtype=odml.DType.float,
                                                 unit='degree'),
                                definition='Rotation of apparatus from '
                                           'vertical axis'))

    parent.append(odml.Property(name='PosX',
                                value=odml.Value(data=5.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='X position of apparatus'))

    parent.append(odml.Property(name='PosY',
                                value=odml.Value(data=13.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='Y position of apparatus'))

    parent.append(odml.Property(name='PosZ',
                                value=odml.Value(data=14.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='Z position of apparatus'))

    parent.append(odml.Property(name='PosOrigin',
                                value=odml.Value(data='midline, waist-level',
                                                 dtype=odml.DType.string),
                                definition='Origin for xyz position of '
                                           'apparatus'))

    parent.append(odml.Property(name='PullRange',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.float,
                                                  unit='mm')
                                       for data in [0.0, 15.0]],
                                definition='Minimum and maximum pull range '
                                           'of target object'))

    parent.append(odml.Property(name='WeightRange',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.float,
                                                  unit='g')
                                       for data in [0.0, 800.0]],
                                definition='Minimum and maximum weight range '
                                           'of target object'))

    parent.append(odml.Property(name='WeightSwitch',
                                value=odml.Value(data='activation of magnet',
                                                 dtype=odml.DType.string),
                                definition='Mechanism to switch weight '
                                           'configuration'))

    parent.append(odml.Property(name='Sensors',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['FSRSensor',
                                                    'HESensor',
                                                    'LoadForce']],
                                definition='Sensors providing additional '
                                           'analogue signals to monitor '
                                           'task and behavior'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['TargetObject']['FSRSensor']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Function',
                                value=odml.Value(data='measures grip and '
                                                      'pulling load forces by '
                                                      'means of force '
                                                      'sensitive resistances '
                                                      '(FSR)',
                                                 dtype=odml.DType.string),
                                definition='Function of apparatus'))

    parent.append(odml.Property(name='SGChannelIDs',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.int)
                                       for data in [138, 140]],
                                definition='Channel IDs of FSR sensors to '
                                           'measure load forces of side grip'))

    parent.append(odml.Property(name='PGChannelIDs',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.int)
                                       for data in [137, 139]],
                                definition='Channel IDs of FSR sensors to '
                                           'measure load forces of precision '
                                           'grip'))

    parent.append(odml.Property(name='SignalType',
                                value=odml.Value(data='analog',
                                                 dtype=odml.DType.string),
                                definition='Signal type of apparatus (digital '
                                           'or analog)'))

    parent.append(odml.Property(name='Events',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['OT', 'OR', 'HS', 'BTB']],
                                definition='Names of events extracted from '
                                           'this signal'))

    parent.append(odml.Property(name='SamplingRate',
                                value=odml.Value(data=1000,
                                                 dtype=odml.DType.int,
                                                 unit='samples/sec'),
                                definition='Sampling rate of data'))

    parent.append(odml.Property(name='ConnectedTo',
                                value=odml.Value(data='NeuralSignalProcessor',
                                                 dtype=odml.DType.string),
                                definition='Target device of apparatus'))

    parent.append(odml.Property(name='MonitoredBy',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Control software of apparatus'))

    parent.append(odml.Property(name='ControlledBy',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Control software of apparatus'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='ns2',
                                                 dtype=odml.DType.string),
                                definition='File in which apparatus signals '
                                           'are saved'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['TargetObject']['HESensor']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Function',
                                value=odml.Value(data='hall effect (HE) sensor'
                                                      ' measures horizontal '
                                                      'displacement',
                                                 dtype=odml.DType.string),
                                definition='Function of apparatus'))

    parent.append(odml.Property(name='Count',
                                value=odml.Value(data=1,
                                                 dtype=odml.DType.int),
                                definition='Number of hall-effect sensors'))

    parent.append(odml.Property(name='ChannelID',
                                value=odml.Value(data=143,
                                                 dtype=odml.DType.int),
                                definition='Channel ID of HESensor to '
                                           'measure horizontal object '
                                           'displacement'))

    parent.append(odml.Property(name='Events',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['OT', 'OR', 'HS', 'BTB']],
                                definition='Names of events extracted from '
                                           'this signal'))

    parent.append(odml.Property(name='SignalType',
                                value=odml.Value(data='analog',
                                                 dtype=odml.DType.string),
                                definition='Signal type of apparatus (digital '
                                           'or analog)'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='ns2',
                                                 dtype=odml.DType.string),
                                definition='File in which apparatus signals '
                                           'are saved'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['TargetObject']['LoadForce']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Function',
                                value=odml.Value(data='used to detect load '
                                                      'force of the target '
                                                      'object',
                                                 dtype=odml.DType.string),
                                definition='Function of apparatus'))

    parent.append(odml.Property(name='Count',
                                value=odml.Value(data=1,
                                                 dtype=odml.DType.int),
                                definition='Number of load force signals'))

    parent.append(odml.Property(name='ChannelID',
                                value=odml.Value(data=141,
                                                 dtype=odml.DType.int),
                                definition='Channel ID of load force signal'))

    parent.append(odml.Property(name='SignalType',
                                value=odml.Value(data='analog',
                                                 dtype=odml.DType.string),
                                definition='Signal type of apparatus (digital '
                                           'or analog)'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='ns2',
                                                 dtype=odml.DType.string),
                                definition='File in which apparatus signals '
                                           'are saved'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['CueSystem']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='CueCodes',
                               type='setup/daq',
                               definition='Information on cue codes'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Function',
                                value=odml.Value(data='instructing movement '
                                                      'configurations of each '
                                                      'trial',
                                                 dtype=odml.DType.string),
                                definition='Function of apparatus'))

    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='visual, LED cube',
                                                 dtype=odml.DType.string),
                                definition='Type of apparatus'))

    parent.append(odml.Property(name='SignalType',
                                value=odml.Value(data='digital',
                                                 dtype=odml.DType.string),
                                definition='Signal type of apparatus (digital '
                                           'or analog)'))

    parent.append(odml.Property(name='EventNames',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['FP-ON', 'CUE-ON',
                                                    'CUE-OFF', 'GO-ON',
                                                    'GO-OFF', 'ERROR-FLASH',
                                                    'VALIDATION-FLASH']],
                                definition='Names of events extracted from '
                                           'this signal'))

    parent.append(odml.Property(name='Length',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Length of apparatus'))

    parent.append(odml.Property(name='Width',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Width of apparatus'))

    parent.append(odml.Property(name='Height',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Height of apparatus'))

    parent.append(odml.Property(name='PosX',
                                value=odml.Value(data=5.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='X position of apparatus'))

    parent.append(odml.Property(name='PosY',
                                value=odml.Value(data=13.0,
                                                 dtype=odml.DType.float,
                                                 unit='cm'),
                                definition='Y position of apparatus'))

    parent.append(odml.Property(name='PosZ',
                                value=default_value(odml.DType.float,
                                                    unit='cm'),
                                definition='Z position of apparatus'))

    parent.append(odml.Property(name='PosOrigin',
                                value=odml.Value(data='midline, waist-level',
                                                 dtype=odml.DType.string),
                                definition='Origin for xyz position of '
                                           'apparatus'))

    parent.append(odml.Property(name='LEDCount',
                                value=odml.Value(data=5,
                                                 dtype=odml.DType.int),
                                definition='Number of LEDs'))

    parent.append(odml.Property(name='LEDConfiguration',
                                value=odml.Value(data='2-1-2',
                                                 dtype=odml.DType.string),
                                definition='Arrangement of LEDs'))

    parent.append(odml.Property(name='CornerLEDColor',
                                value=odml.Value(data='red',
                                                 dtype=odml.DType.string),
                                definition='Color of the corner LEDs'))

    parent.append(odml.Property(name='CenterLEDColor',
                                value=odml.Value(data='yellow',
                                                 dtype=odml.DType.string),
                                definition='Color of the center LEDs'))

    parent.append(odml.Property(name='SamplingRate',
                                value=odml.Value(data=30000,
                                                 dtype=odml.DType.int,
                                                 unit='samples/sec'),
                                definition='Sampling rate of data'))

    parent.append(odml.Property(name='ConnectedTo',
                                value=odml.Value(data='NeuralSignalProcessor',
                                                 dtype=odml.DType.string),
                                definition='Target device of apparatus'))

    parent.append(odml.Property(name='MonitoredBy',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Control software of apparatus'))

    parent.append(odml.Property(name='ControlledBy',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Control software of apparatus'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='nev',
                                                 dtype=odml.DType.string),
                                definition='File in which apparatus signals '
                                           'are saved'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['CueSystem']['CueCodes']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='LeftLEDs',
                                value=odml.Value(data='side grip',
                                                 dtype=odml.DType.string),
                                definition='LEDs coding this trial '
                                           'instruction'))

    parent.append(odml.Property(name='RightLEDs',
                                value=odml.Value(data='precision grip',
                                                 dtype=odml.DType.string),
                                definition='LEDs coding this trial '
                                           'instruction'))

    parent.append(odml.Property(name='BottomLEDs',
                                value=odml.Value(data='low force',
                                                 dtype=odml.DType.string),
                                definition='LEDs coding this trial '
                                           'instruction'))

    parent.append(odml.Property(name='TopLEDs',
                                value=odml.Value(data='high force',
                                                 dtype=odml.DType.string),
                                definition='LEDs coding this trial '
                                           'instruction'))

    parent.append(odml.Property(name='CenterLED',
                                value=odml.Value(data='fixation point',
                                                 dtype=odml.DType.string),
                                definition='LEDs coding this trial '
                                           'instruction'))

    parent.append(odml.Property(name='CornerLEDs',
                                value=odml.Value(data='grip error',
                                                 dtype=odml.DType.string),
                                definition='LEDs coding this trial '
                                           'instruction'))

    # SET PARENT NODE
    parent = doc['Setup']['Apparatus']['RewardSystem']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Function',
                                value=odml.Value(data='rewards subject if '
                                                      'trial was performed '
                                                      'correctly',
                                                 dtype=odml.DType.string),
                                definition='Function of apparatus'))

    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='reward pipe',
                                                 dtype=odml.DType.string),
                                definition='Type of apparatus'))

    parent.append(odml.Property(name='SignalType',
                                value=odml.Value(data='digital',
                                                 dtype=odml.DType.string),
                                definition='Signal type of apparatus (digital '
                                           'or analog)'))

    parent.append(odml.Property(name='EventNames',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['RW']],
                                definition='Names of events extracted from '
                                           'this signal'))

    parent.append(odml.Property(name='RewardFood',
                                value=odml.Value(data='apple sauce',
                                                 dtype=odml.DType.string),
                                definition='Food used for reward'))

    parent.append(odml.Property(name='RewardAmount',
                                value=odml.Value(data=0.5,
                                                 dtype=odml.DType.float,
                                                 unit='ml/200ms'),
                                definition='Food amount used for reward per '
                                           'trial'))

    parent.append(odml.Property(name='SamplingRate',
                                value=odml.Value(data=30000,
                                                 dtype=odml.DType.int,
                                                 unit='samples/sec'),
                                definition='Sampling rate of data'))

    parent.append(odml.Property(name='ConnectedTo',
                                value=odml.Value(data='NeuralSignalProcessor',
                                                 dtype=odml.DType.string),
                                definition='Target device of apparatus'))

    parent.append(odml.Property(name='MonitoredBy',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Control software of apparatus'))

    parent.append(odml.Property(name='ControlledBy',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Control software of apparatus'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='nev',
                                                 dtype=odml.DType.string),
                                definition='File in which apparatus signals '
                                           'are saved'))

    # SET PARENT NODE
    parent = doc['Setup']['ControlComputer']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='ControlSoftware',
                               type='setup/daq',
                               definition='Information on the control software'
                                          ' of the Cerebus^TM system'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='CPU',
                                value=default_value(odml.DType.string),
                                definition='CPU model'))

    parent.append(odml.Property(name='RAM',
                                value=default_value(odml.DType.float,
                                                    unit='GB'),
                                definition='Amount of RAM'))

    parent.append(odml.Property(name='Mainboard',
                                value=default_value(odml.DType.string),
                                definition='Mainboard type'))

    parent.append(odml.Property(name='HarddiskCapacity',
                                value=default_value(odml.DType.float,
                                                    unit='GB'),
                                definition='Disk capacity of data storage'))

    parent.append(odml.Property(name='OSType',
                                value=default_value(odml.DType.string),
                                definition='Type of operating system'))

    parent.append(odml.Property(name='OSVersion',
                                value=default_value(odml.DType.string),
                                definition='Version of operating system'))

    # SET PARENT NODE
    parent = doc['Setup']['ControlComputer']['ControlSoftware']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Name',
                                value=odml.Value(data='LabVIEW',
                                                 dtype=odml.DType.string),
                                definition='Name of software'))

    parent.append(odml.Property(name='Manufacturer',
                                value=odml.Value(data='National Instruments',
                                                 dtype=odml.DType.string),
                                definition='Manufacturer of software'))

    parent.append(odml.Property(name='Version',
                                value=default_value(odml.DType.string),
                                definition='Software version'))

    return doc
