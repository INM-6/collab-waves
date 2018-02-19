# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 16:30:05 2014

@author: zehl
"""
import odml
from rgodml.odml_io import default_value


def get_DIO_subsection(DIO_port):
    """
    Constructs a Section describing the specified DIO_port including
    Properties, Values and (possibly default) data

    Args:
        DIO_port (string):
            Type of the digital input/output port

    Returns:
        (odml.section.BaseSection)
            Section describing the given DIO_port

    Note:
        Valid digital input/output ports: ExpDigitalInput, ExpDigitalOutput,
        SerialDigitalIO, ExtSync, NSPSync
    """
    if DIO_port == 'ExpDigitalInput':
        # SET PARENT NODE
        parent = odml.Section(name='ExpDigitalInput',
                              type='setup/daq/hardware',
                              definition='Information on the experimental '
                                         'digital inputs')

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='PortType',
                                    value=odml.Value(data='DB-37',
                                                     dtype=odml.DType.string),
                                    definition='Connector port type'))

        parent.append(odml.Property(name='ChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Total number of input '
                                               'channels'))

        parent.append(odml.Property(name='UsedChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Number of used input channels '
                                               'usded'))

        parent.append(odml.Property(name='DataBits',
                                    value=odml.Value(data=16,
                                                     dtype=odml.DType.int,
                                                     unit='bit'),
                                    definition='Width of digital data in '
                                               'bits'))

        parent.append(odml.Property(name='SamplingRate',
                                    value=odml.Value(data=30000,
                                                     dtype=odml.DType.int,
                                                     unit='samples/sec'),
                                    definition='Sampling rate of data'))

        parent.append(odml.Property(name='InFrom',
                                    value=odml.Value(data='-',
                                                     dtype=odml.DType.string),
                                    definition='Devices from which inputs are '
                                               'received'))

        parent.append(odml.Property(name='SavedIn',
                                    value=odml.Value(data='nev',
                                                     dtype=odml.DType.string),
                                    definition='File in which data are saved'))

    elif DIO_port == 'ExpDigitalOutput':
        # SET PARENT NODE
        parent = odml.Section(name='ExpDigitalOutput',
                              type='setup/daq/hardware',
                              definition='Information on the experimental '
                                         'digital outputs')

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='PortType',
                                    value=odml.Value(data='BNC',
                                                     dtype=odml.DType.string),
                                    definition='Connector port type'))

        parent.append(odml.Property(name='ChannelCount',
                                    value=odml.Value(data=4,
                                                     dtype=odml.DType.int),
                                    definition='Total number of output '
                                               'channels'))

        parent.append(odml.Property(name='UsedChannelCount',
                                    value=default_value(odml.DType.int),
                                    definition='Number of output channels '
                                               'used'))

        parent.append(odml.Property(name='DataBits',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int,
                                                     unit='bit'),
                                    definition='Width of digital data in '
                                               'bits'))

        parent.append(odml.Property(name='Frequency',
                                    value=default_value(odml.DType.float,
                                                        unit='Hz'),
                                    definition='Frequency of the user-defined '
                                               'pulse train'))

        parent.append(odml.Property(name='DutyCycle',
                                    value=default_value(odml.DType.float,
                                                        unit='%'),
                                    definition='Ratio of pulse duration to '
                                               'period of the user-defined '
                                               'pulse train'))

        parent.append(odml.Property(name='OutTo',
                                    value=default_value(odml.DType.string),
                                    definition='Target device of each output'))

        parent.append(odml.Property(name='SavedIn',
                                    value=default_value(odml.DType.string),
                                    definition='File in which data are saved'))

    elif DIO_port == 'SerialDigitalInput':
        # SET PARENT NODE
        parent = odml.Section(name='SerialDigitalInput',
                              type='setup/daq/hardware',
                              definition='Information on the experimental '
                                         'serial digital input')

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='PortType',
                                    value=odml.Value(data='RS232 DB-9',
                                                     dtype=odml.DType.string),
                                    definition='Connector port type'))

        parent.append(odml.Property(name='ChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Total number of input '
                                               'channels'))

        parent.append(odml.Property(name='UsedChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Number of input channels '
                                               'used'))

        parent.append(odml.Property(name='DataBits',
                                    value=odml.Value(data=8,
                                                     dtype=odml.DType.int,
                                                     unit='bit'),
                                    definition='Width of digital data in '
                                               'bits'))

        parent.append(odml.Property(name='BaudRate',
                                    value=odml.Value(data=115200,
                                                     dtype=odml.DType.int,
                                                     unit='bit/sec'),
                                    definition='Baud rate of digital input'))

        parent.append(odml.Property(name='Parity',
                                    value=odml.Value(data=False,
                                                     dtype=odml.DType.boolean),
                                    definition='Data parity bit set'))

        parent.append(odml.Property(name='StopBits',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int,
                                                     unit='bit'),
                                    definition='Number of stop bits'))

        parent.append(odml.Property(name='FlowControl',
                                    value=odml.Value(data='disabled',
                                                     dtype=odml.DType.string),
                                    definition='Flow control specification'))

        parent.append(odml.Property(name='InFrom',
                                    value=default_value(odml.DType.string),
                                    definition='Devices from which inputs are '
                                               'received'))

        parent.append(odml.Property(name='SavedIn',
                                    value=default_value(odml.DType.string),
                                    definition='File in which data are saved'))

    elif DIO_port == 'SerialDigitalOutput':
        # SET PARENT NODE
        parent = odml.Section(name='SerialDigitalOutput',
                              type='setup/daq/hardware',
                              definition='Information on this experimental '
                                         'serial digital output')

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='PortType',
                                    value=odml.Value(data='RS232 DB-9',
                                                     dtype=odml.DType.string),
                                    definition='Connector port type'))

        parent.append(odml.Property(name='ChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Total number of output '
                                               'channels'))

        parent.append(odml.Property(name='UsedChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Number of output channels '
                                               'used'))

        parent.append(odml.Property(name='DataBits',
                                    value=odml.Value(data=8,
                                                     dtype=odml.DType.int,
                                                     unit='bit'),
                                    definition='Width of digital data in '
                                               'bits'))

        parent.append(odml.Property(name='BaudRate',
                                    value=odml.Value(data=115200,
                                                     dtype=odml.DType.int,
                                                     unit='bit/sec'),
                                    definition='Baud rate of digital input'))

        parent.append(odml.Property(name='Parity',
                                    value=odml.Value(data=False,
                                                     dtype=odml.DType.boolean),
                                    definition='Data parity bit set'))

        parent.append(odml.Property(name='StopBits',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int,
                                                     unit='bit'),
                                    definition='Number of stop bits'))

        parent.append(odml.Property(name='FlowControl',
                                    value=odml.Value(data='disabled',
                                                     dtype=odml.DType.string),
                                    definition='Flow control specification'))

        parent.append(odml.Property(name='OutTo',
                                    value=default_value(odml.DType.string),
                                    definition='Target device of each output'))

        parent.append(odml.Property(name='SavedIn',
                                    value=default_value(odml.DType.string),
                                    definition='File in which data are saved'))

    elif DIO_port == 'ExtSync':
        # SET PARENT NODE
        parent = odml.Section(name='ExtSync',
                              type='setup/daq/hardware',
                              definition='Information on the synchronization '
                                         'port used to inform external devices'
                                         ' when the neural signal inputs and '
                                         'front panel ports are scanned')

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='PortType',
                                    value=odml.Value(data='BCN',
                                                     dtype=odml.DType.string),
                                    definition='Connector port type'))

        parent.append(odml.Property(name='ChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Total number of ports of this '
                                               'type'))

        parent.append(odml.Property(name='UsedChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Number of used ports of this '
                                               'type'))

        parent.append(odml.Property(name='Trigger',
                                    value=odml.Value(data='rising edge',
                                                     dtype=odml.DType.string),
                                    definition='Activation condition'))

        parent.append(odml.Property(name='OutTo',
                                    value=default_value(odml.DType.string),
                                    definition='Target device of each output'))

    elif DIO_port == 'NSPSync':
        # SET PARENT NODE
        parent = odml.Section(name='NSPSync',
                              type='setup/daq/hardware',
                              definition='Information on the port used to '
                                         'synchronize two NSPs')

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='PortType',
                                    value=odml.Value(data='DB-9',
                                                     dtype=odml.DType.string),
                                    definition='Connector port type'))

        parent.append(odml.Property(name='ChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Total number of ports of this '
                                               'type'))

        parent.append(odml.Property(name='UsedChannelCount',
                                    value=odml.Value(data=1,
                                                     dtype=odml.DType.int),
                                    definition='Number of used ports of this '
                                               'type'))

        parent.append(odml.Property(name='PairedNSP',
                                    value=default_value(odml.DType.string),
                                    definition='ID of NSP synchronized with'))

    else:
        raise TypeError('Unknown DIO_port')

    return parent


def get_connector_subsection(connector_style=None):
    """
    Constructs a Section describing the a connector of the specified
    connector_style including Properties, Values and (possible default) data

    If the connector_style is None, it returns a default Section for a
    connector of any connector_style

    Args:
        connector_style (string, optional):
            Connector style of the used connector

    Returns:
        (odml.section.BaseSection):
            Section describing a connector of the given connector_style

    Note:
        Valid connector_styles: CerePort, ICS-96, Omnetics or None (default)
    """
    # SET PARENT NODE
    parent = odml.Section(name='Connector',
                          type='setup/daq/hardware',
                          definition='Information on the used connector')

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Style',
                                value=odml.Value(data=connector_style,
                                                 dtype=odml.DType.string)
                                if connector_style
                                else default_value(odml.DType.string),
                                definition='Style of the connector'))

    parent.append(odml.Property(name='Height',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Height of the connector'))

    parent.append(odml.Property(name='DiameterBase',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Diameter of the base of the '
                                           'connector'))

    parent.append(odml.Property(name='DiameterBody',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Diameter of the body of the '
                                           'connector'))

    parent.append(odml.Property(name='Width',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Width of the connector'))

    parent.append(odml.Property(name='Length',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Length of the connector'))

    parent.append(odml.Property(name='ImplantType',
                                value=default_value(odml.DType.string),
                                definition='Type of implantation'))

    parent.append(odml.Property(name='SkullMount',
                                value=default_value(odml.DType.string),
                                definition='Type of mounting to the skull'))

    if connector_style == 'CerePort':
        parent.properties['Height'].value.data = 16.5
        parent.properties['DiameterBase'].value.data = 19.0
        parent.properties['DiameterBody'].value.data = 11.0
        parent.properties['ImplantType'].value.data = 'chronic'
        parent.properties['SkullMount'].value.data = 'bone screws'
        parent.remove(parent.properties['Width'])
        parent.remove(parent.properties['Length'])

    elif connector_style == 'ICS-96':
        parent.properties['Height'].value.data = 18.0
        parent.properties['Width'].value.data = 17.0
        parent.properties['Length'].value.data = 37.0
        parent.properties['ImplantType'].value.data = 'acute'
        parent.properties['SkullMount'].value.data = 'NA'
        parent.remove(parent.properties['DiameterBase'])
        parent.remove(parent.properties['DiameterBody'])

    elif connector_style == 'Omnetics':
        parent.properties['Height']['Val'].data = 9.0
        parent.properties['Width']['Val'].data = 7.0
        parent.properties['Length']['Val'].data = 13.0
        parent.properties['ImplantType']['Val'].data = 'chronic'
        parent.properties['SkullMount']['Val'].data = 'dental acrylic'
        parent.remove(parent.properties['DiameterBase'])
        parent.remove(parent.properties['DiameterBody'])

    elif connector_style is None:
        pass

    else:
        raise TypeError('Unknown connector_style')

    return parent


def get_grid_section_names(grid_no=0):
    """
    Returns a list of grid Section names

    Args:
        grid_no (int, optional):
            Number of grids

    Returns:
        (list of strings):
            List of grid Section names

    Note:
        Grid Section names are formatted as 'Grid_ID' where ID is a given as a
        two digit list index (e.g. ['Grid_01','Grid_02'] for grid_no = 2)

        If grid_no = 0 the list ['Grid_XX'] is returned
    """
    if grid_no == 0:
        section_names = ['Grid_XX']

    elif grid_no > 0:
        section_names = ['Grid_%02i' % i for i in range(1, grid_no+1)]

    else:
        raise ValueError('Grid number must be non-negative')

    return section_names


def get_elec_section_names(electrode_ids=[]):
    """
    Returns a list of electrode Section names

    Args:
        electrode_no (int, optional):
            Number of electrodes

    Returns:
        (list of strings):
            List of electrode Section names

    Note:
        Electrode Section names are formatted as 'Electrode_ID' where ID is a
        given as a three digit list index (e.g. ['Electrode_001',
        'Electrode_002'] for electrode_no = 2)

        If electrode_no = 0 the list ['Electrode_XXX'] is returned
    """
    if len(electrode_ids) == 0:
        section_names = ['Electrode_XXX']

    elif len(electrode_ids) > 0:
        section_names = ['Electrode_%03i' % i for i in electrode_ids]

    else:
        raise ValueError('Electrode number must be non-negative')

    return section_names


def get_NSIO_filter_sets(avail_nsx):
    """
    Returns a dict with filter Section names as keys and the corresponding
    filter settings as subkeys

    Args:
        avail_nsx (list of int or str):
            Number of electrodes

    Returns:
        (dict):
            Dictionary with filter Section names as keys, and corresponding
            filter settings as subkeys. Possible are filters for ns1 to ns6
            with the following settings (as odml.Values):
                nsx:
                    HighPassFreq:
                        default_value(odml.DType.float, unit='Hz')
                    LowPassFreq:
                        default_value(odml.DType.float, unit='Hz')
                    HighPassOrder:
                        default_value(odml.DType.int)
                    LowPassOrder:
                        default_value(odml.DType.int)
                    SavedIn:
                        odml.Value(nsx, odml.DType.string)
    Note:
        The nsx is defined as 'Filter_nsX' where X is given as a digit (e.g.
        ['Filter_ns2'] for avail_nsx = [2])

        If avail_nsx = [] the settings for 'Filter_nsX' are returned
    """
    pfsets = {'nsX': {'hpfr': -1.0, 'lpfr': -1.0, 'hpor': -1, 'lpor': -1},
              'ns2': {'hpfr': 0.3, 'lpfr': 250., 'hpor': -1, 'lpor': -1},
              'ns5': {'hpfr': 500., 'lpfr': 7500., 'hpor': -1, 'lpor': -1},
              'ns6': {'hpfr': -1.0, 'lpfr': -1.0, 'hpor': -1, 'lpor': -1}}

    fsets = {}
    if len(avail_nsx) == 0:
        fsets['Filter_nsX'] = {}

    elif len(avail_nsx) > 0:
        if isinstance(avail_nsx[0], str):
            fsets = dict(('Filter_%s' % i, {}) for i in avail_nsx)
        elif isinstance(avail_nsx[0], int):
            fsets = dict(('Filter_ns%i' % i, {}) for i in avail_nsx)
        else:
            raise ValueError('avail_nsX must be a list of str or int')

    for f in fsets.keys():
        nsx = f[-3:]
        fsets[f]['HighPassFreq'] = odml.Value(data=pfsets[nsx]['hpfr'],
                                              dtype=odml.DType.float,
                                              unit='Hz')
        fsets[f]['LowPassFreq'] = odml.Value(data=pfsets[nsx]['lpfr'],
                                             dtype=odml.DType.float,
                                             unit='Hz')
        fsets[f]['HighPassOrder'] = odml.Value(data=pfsets[nsx]['hpor'],
                                               dtype=odml.DType.int)
        fsets[f]['LowPassOrder'] = odml.Value(data=pfsets[nsx]['lpor'],
                                              dtype=odml.DType.int)

    return fsets


def get_headstage_sets(used_headstage):
    """
    Returns a dict with the given settings for the headstage Section

    Args:
        used_headstage (str):
            Model name of used headstage

    Returns:
        (dict):
            Dictionary with headstage settings (default settings):
                Model:
                    default_value(odml.DType.string)
                Type:
                    default_value(odml.DType.string)
                Gain:
                    default_value(odml.DType.int)
                InFrom:
                    default_value(odml.DType.string)
                OutTo:
                    default_value(odml.DType.string)
    Note:
        Possible headstage models: Samtec, Omnetics, PatientCable, CerePlexE,
        CerePlexM

        If used_headstage = None default settings are returned
    """
    hmodels = {}
    for m in ['Samtec', 'Omnetics', 'PatientCable']:
        hmodels[m] = {'Type': 'analog', 'InFrom': 'UtahArray',
                      'OutTo': 'NeuralSignalAmplifier'}
    for m in ['CerePlexE', 'CerePlexM']:
        hmodels[m] = {'Type': 'digital', 'InFrom': 'UtahArray',
                      'OutTo': 'DigitalHub'}
    hmodels['-'] = {'Type': '-', 'InFrom': '-', 'OutTo': '-'}

    if not used_headstage:
        used_headstage = '-'

    hsets = {'Model': odml.Value(data=used_headstage,
                                 dtype=odml.DType.string),
             'Type': odml.Value(data=hmodels[used_headstage]['Type'],
                                dtype=odml.DType.string),
             'InFrom': odml.Value(data=hmodels[used_headstage]['InFrom'],
                                  dtype=odml.DType.string),
             'OutTo': odml.Value(data=hmodels[used_headstage]['OutTo'],
                                 dtype=odml.DType.string)}

    return hsets


def blackrock_tree(used_DIOports=[], connector_style=None, grid_no=0,
                   electrode_ids=[], avail_nsx=[], analogI_chIDs=[],
                   used_headstage=None, neural_signal_stimulator=False,
                   spike_sorted=False):
    """
    Builds an odML Document for Blackrock metadata

    Args:
        used_DIOports (list of strings, optional):
            List of digital input and output ports used in the setup
        connector_style (string, optional):
            Style of connector used in the setup (If not specified a default
            connector Section is generated)
        grid_no (int):
            Number of array grids connected to the specified connector (If
            grid_no = 0 a single default grid Section is generated)
        electrode_ids (list of int):
            IDs of active electrodes (across all grids) (If electrode_ids = []
            a single default electrode Section is generated)
        avail_nsx (list of int or str):
            Available nsX files (given either as 'nsX' or as int) (If
            avail_nsx = [] a default Sections and Values are generated)
        analogI_chIDs (list of int):
            List of channel IDs of analog inputs recorded and saved in ns2 (If
            analogI_chIDs = [] default analog input Values are generated)
        used_headstage (str):
            Name of used headstage
    Returns:
        (odml.doc.BaseDocument)
            odML Document with default Blackrock metadata
    """

    # CREATE A DOCUMENT
    doc = odml.Document(version='0.1')

    # CREATE A MAIN SECTION
    doc.append(odml.Section(name='Cerebus',
                            type='setup/daq',
                            definition='Information on the Cerebus^TM data '
                                       'acquisition system'))

    doc.append(odml.Section(name='UtahArray',
                            type='setup/daq/hardware',
                            definition='Information on the Blackrock Utah '
                                       'Electrode Array'))

    doc.append(odml.Section(name='Headstage',
                            type='setup/daq/hardware',
                            definition='Information on the headstage used'))

    # SET PARENT NODE
    parent = doc['Cerebus']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='NeuralSignalProcessor',
                               type='setup/daq/hardware',
                               definition='Information of the Cerebus^TM '
                                          'neural signal processor (NSP)'))

    parent.append(odml.Section(name='NeuralSignalAmplifier',
                               type='setup/daq/hardware',
                               definition='Information on the Cerebus^TM '
                                          'neural signal amplifier (NSA)'))

    if neural_signal_stimulator:
        parent.append(odml.Section(name='NeuralSignalStimulator',
                                   type='setup/daq/hardware',
                                   definition='Information on the Cerebus^TM '
                                              'neural signal stimulator '
                                              '(NSS)'))

    parent.append(odml.Section(name='ControlComputer',
                               type='setup/daq/software',
                               definition='Information on the Cerebus^TM '
                                          'control computer'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Owner',
                                value=odml.Value(data='Riken, Japan',
                                                 dtype=odml.DType.string),
                                definition='Owner of the device'))

    parent.append(odml.Property(name='Manufacturer',
                                value='Blackrock Microsystems',
                                definition='Manufacturer of the device'))

    parent.append(odml.Property(name='UserManual',
                                value=default_value(odml.DType.url),
                                definition='User manual URL of the device'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalProcessor']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='AnalogSignals',
                               type='setup/daq',
                               definition='Information on external analog '
                                          'signals from or to the device'))

    parent.append(odml.Section(name='NeuralSignals',
                               type='setup/daq',
                               definition='Information on neural signals '
                                          'of the device'))

    parent.append(odml.Section(name='DigitalSignals',
                               type='setup/daq',
                               definition='Information on digital signals '
                                          'from and to the device'))

    parent.append(odml.Section(name='SpikeDetection',
                               type='setup/daq',
                               definition='Information on spike detection'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='Real-time',
                                                 dtype=odml.DType.string),
                                definition='Type of the device'))

    parent.append(odml.Property(name='Model',
                                value=default_value(odml.DType.string),
                                definition='Model of the device'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalProcessor']['AnalogSignals']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='ADConverter',
                               type='setup/daq',
                               definition='Information on the analog-to-'
                                          'digital converter'))

    parent.append(odml.Section(name='Input',
                               type='setup/daq',
                               definition='Information on analog input '
                                          'signals'))

    parent.append(odml.Section(name='Output',
                               type='setup/daq',
                               definition='Information on analog output '
                                          'signals'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='UsedPortType',
                                value=odml.Value(data='BNC',
                                                 dtype=odml.DType.string),
                                definition='Connector port type'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalProcessor']['AnalogSignals']['Input']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='ChannelCount',
                                value=odml.Value(data=16,
                                                 dtype=odml.DType.int),
                                definition='Total number of input channels'))

    parent.append(odml.Property(name='UsedChannelCount',
                                value=odml.Value(data=len(analogI_chIDs),
                                                 dtype=odml.DType.int)
                                if len(analogI_chIDs) > 0
                                else default_value(odml.DType.int),
                                definition='Number of input channels used'))

    parent.append(odml.Property(name='ACChannelCount',
                                value=odml.Value(data=8,
                                                 dtype=odml.DType.int),
                                definition='Number of AC coupled input '
                                           'channels (IDs 1 - 9)'))

    parent.append(odml.Property(name='UsedACChannelCount',
                                value=odml.Value(data=0,
                                                 dtype=odml.DType.int),
                                definition='Number of used AC coupled input '
                                           'channels'))

    parent.append(odml.Property(name='DCChannelCount',
                                value=odml.Value(data=8,
                                                 dtype=odml.DType.int),
                                definition='Number of DC coupled input '
                                           'channels (IDs 9 -16)'))

    parent.append(odml.Property(name='UsedDCChannelCount',
                                value=odml.Value(data=8,
                                                 dtype=odml.DType.int),
                                definition='Number of used DC coupled input '
                                           'channels'))

    parent.append(odml.Property(name='InFrom',
                                value=odml.Value(data='TargetObject',
                                                 dtype=odml.DType.string),
                                definition='Devices from which inputs are '
                                           'received'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='ns2',
                                                 dtype=odml.DType.string),
                                definition='File in which input data are '
                                           'saved'))

    parent.append(odml.Property(name='SavedChannelIDs',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.int)
                                       for data in analogI_chIDs]
                                if len(analogI_chIDs) > 0
                                else default_value(odml.DType.int),
                                definition='File in which input data are '
                                           'saved'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalProcessor']['AnalogSignals']['Output']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='ChannelCount',
                                value=odml.Value(data=4,
                                                 dtype=odml.DType.int),
                                definition='Total number of output channels'))

    parent.append(odml.Property(name='UsedChannelCount',
                                value=odml.Value(data=0,
                                                 dtype=odml.DType.int),
                                definition='Number of output channels used'))

    parent.append(odml.Property(name='OutTo',
                                value=default_value(odml.DType.string),
                                definition='Target device of each output'))

    parent.append(odml.Property(name='SavedIn',
                                value=default_value(odml.DType.string),
                                definition='File in which output data are '
                                           'saved'))

    # SET PARENT NODE
    parent = \
        doc['Cerebus']['NeuralSignalProcessor']['AnalogSignals']['ADConverter']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='AIRange',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.float,
                                                  unit='V')
                                       for data in [-5.0, 5.0]],
                                definition='Voltage range of analog input'))

    parent.append(odml.Property(name='MaxAIImpedance',
                                value=odml.Value(data=100.0,
                                                 dtype=odml.DType.float,
                                                 unit=u"\u2126"),
                                definition='Maximum impedance of analog '
                                           'input'))

    parent.append(odml.Property(name='Resolution',
                                value=odml.Value(data=0.15,
                                                 dtype=odml.DType.float,
                                                 unit='mV/bit'),
                                definition='Resolution of AD converter'))

    parent.append(odml.Property(name='DataBits',
                                value=odml.Value(data=16,
                                                 dtype=odml.DType.int,
                                                 unit='bit'),
                                definition='Width of digital data in bits'))

    parent.append(odml.Property(name='SamplingRate',
                                value=odml.Value(data=1000,
                                                 dtype=odml.DType.int,
                                                 unit='samples/sec'),
                                definition='Sampling rate of data'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalProcessor']['NeuralSignals']

    # APPEND SUBSECTIONS
    for sname in get_NSIO_filter_sets(avail_nsx).keys():
        parent.append(odml.Section(name=sname,
                                   type='setup/daq',
                                   definition='Information on the filter '
                                              'settings'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='UsedPortType',
                                value=odml.Value(data='Fiber-optic link',
                                                 dtype=odml.DType.string),
                                definition='Digital port type'))

    parent.append(odml.Property(name='ChannelCount',
                                value=odml.Value(data=128,
                                                 dtype=odml.DType.int),
                                definition='Total number of input '
                                           'channels'))

    parent.append(odml.Property(name='UsedChannelCount',
                                value=odml.Value(data=len(electrode_ids),
                                                 dtype=odml.DType.int)
                                if len(electrode_ids) > 0
                                else default_value(odml.DType.int),
                                definition='Number of input channels '
                                           'used'))

    parent.append(odml.Property(name='InFrom',
                                value=odml.Value(data='NeuralSignalAmplifier',
                                                 dtype=odml.DType.string),
                                definition='Devices from which inputs are '
                                           'received'))

    parent.append(odml.Property(name='SavedIn',
                                value=[odml.Value(data=data[-3:],
                                                  dtype=odml.DType.string)
                                       for data in
                                       get_NSIO_filter_sets(avail_nsx).keys()],
                                definition='File in which input data are '
                                           'saved'))

    parent.append(odml.Property(name='SavedChannelIDs',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.int)
                                       for data in electrode_ids]
                                if len(electrode_ids) > 0
                                else default_value(odml.DType.int),
                                definition='File in which input data are '
                                           'saved'))

    # SET PARENT NODE
    for sname, fsets in get_NSIO_filter_sets(avail_nsx).items():
        parent = \
            doc['Cerebus']['NeuralSignalProcessor']['NeuralSignals'][sname]

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Type',
                                    value=odml.Value(data='butterworth',
                                                     dtype=odml.DType.string),
                                    definition='Used filter type'))

        parent.append(odml.Property(name='Causal',
                                    value=odml.Value(data=True,
                                                     dtype=odml.DType.boolean),
                                    definition='States if the used filter is '
                                               'causal'))

        parent.append(odml.Property(name='HighPassFreq',
                                    value=fsets['HighPassFreq'],
                                    definition='High pass frequency used'))

        parent.append(odml.Property(name='LowPassFreq',
                                    value=fsets['LowPassFreq'],
                                    definition='Low pass frequency used'))

        parent.append(odml.Property(name='HighPassOrder',
                                    value=fsets['HighPassOrder'],
                                    definition='High pass order used'))

        parent.append(odml.Property(name='LowPassOrder',
                                    value=fsets['LowPassOrder'],
                                    definition='Low pass order used'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalProcessor']['DigitalSignals']

    # APPEND SUBSECTIONS
    for DIO_port in used_DIOports:
        DIOport_sec = get_DIO_subsection(DIO_port)
        parent.append(DIOport_sec)

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Types',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in ['ExpDigitalInput',
                                                    'ExpDigitalOutput',
                                                    'SerialDigitalInput',
                                                    'SerialDigitalOutput',
                                                    'ExtSync', 'NSPSync']],
                                definition='Possible digital input and output '
                                           '(DIO) types'))

    parent.append(odml.Property(name='UsedTypes',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in used_DIOports]
                                if len(used_DIOports) > 0
                                else default_value(odml.DType.string),
                                definition='Used digital input and output '
                                           '(DIO) types'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalProcessor']['SpikeDetection']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='Filter',
                               type='setup/daq',
                               definition='Information on the filter '
                                          'settings'))

    parent.append(odml.Section(name='Waveforms',
                               type='setup/daq',
                               definition='Information on the recorded spike '
                                          'waveforms of the spike detection'))

    parent.append(odml.Section(name='Events',
                               type='setup/daq',
                               definition='Information on the recorded spike '
                                          'events of the spike detection'))

    parent.append(odml.Section(name='Rejections',
                               type='setup/daq',
                               definition='Information on rejected recorded '
                                          'spike events of the spike '
                                          'detection'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='ThresholdType',
                                value=default_value(odml.DType.string),
                                definition='Used method to define a data '
                                           'threshold for the spike '
                                           'detection'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='nev',
                                                 dtype=odml.DType.string),
                                definition='File in which input data are '
                                           'saved'))

    # SET PARENT NODE
    parent = \
        doc['Cerebus']['NeuralSignalProcessor']['SpikeDetection']['Filter']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='butterworth',
                                                 dtype=odml.DType.string),
                                definition='Used filter type'))

    parent.append(odml.Property(name='Causal',
                                value=odml.Value(data=True,
                                                 dtype=odml.DType.boolean),
                                definition='States if the used filter is '
                                           'causal'))

    parent.append(odml.Property(name='HighPassFreq',
                                value=default_value(odml.DType.float,
                                                    unit='Hz'),
                                definition='High pass frequency used'))

    parent.append(odml.Property(name='LowPassFreq',
                                value=default_value(odml.DType.float,
                                                    unit='Hz'),
                                definition='Low pass frequency used'))

    parent.append(odml.Property(name='HighPassOrder',
                                value=default_value(odml.DType.int),
                                definition='High pass order used'))

    parent.append(odml.Property(name='LowPassOrder',
                                value=default_value(odml.DType.int),
                                definition='Low pass order used'))

    # SET PARENT NODE
    parent = \
        doc['Cerebus']['NeuralSignalProcessor']['SpikeDetection']['Waveforms']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Width',
                                value=odml.Value(data=1.6,
                                                 dtype=odml.DType.float,
                                                 unit='ms'),
                                definition='Waveform window width'))

    parent.append(odml.Property(name='SamplingRate',
                                value=odml.Value(data=30000,
                                                 dtype=odml.DType.int,
                                                 unit='samples/sec'),
                                definition='Sampling rate of data'))

    parent.append(odml.Property(name='PreThresholdSamples',
                                value=odml.Value(data=10,
                                                 dtype=odml.DType.int,
                                                 unit='samples'),
                                definition='Number of samples before waveform '
                                           'passes threshold'))

    # SET PARENT NODE
    parent = \
        doc['Cerebus']['NeuralSignalProcessor']['SpikeDetection']['Events']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='SamplingRate',
                                value=odml.Value(data=30000,
                                                 dtype=odml.DType.int,
                                                 unit='samples/sec'),
                                definition='Sampling rate of data'))

    parent.append(odml.Property(name='PreSortedUnitIDRange',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.int)
                                       for data in [1, 16]],
                                definition='Possible unit ID range of online '
                                           'sorted spike events on each '
                                           'channel'))

    parent.append(odml.Property(name='UnclassifiedUnitID',
                                value=odml.Value(data=0,
                                                 dtype=odml.DType.int),
                                definition='Unit ID of unclassified spike '
                                           'events on each channel'))

    # SET PARENT NODE
    parent = \
        doc['Cerebus']['NeuralSignalProcessor']['SpikeDetection']['Rejections']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='AmplitudeLimits',
                                value=default_value(odml.DType.float,
                                                    unit='mV'),
                                definition='Amplitude limits of spike '
                                           'waveforms used to reject spike '
                                           'events'))

    parent.append(odml.Property(name='SyncChannelLimit',
                                value=default_value(odml.DType.int),
                                definition='Minimum number of channels with a '
                                           'synchronously detected spike '
                                           'events used to reject these spike '
                                           'events'))

    parent.append(odml.Property(name='SyncChannelRefPeriod',
                                value=default_value(odml.DType.float,
                                                    unit='samples'),
                                definition='Number of samples for which spike '
                                           'detection is turned off after a '
                                           'rejected synchronous spike event'))

    parent.append(odml.Property(name='UnitID',
                                value=odml.Value(data=255,
                                                 dtype=odml.DType.int),
                                definition='Unit ID reserved for rejected '
                                           'spike events on each channel'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalAmplifier']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='ADConverter',
                               type='setup/daq',
                               definition='Information on the analog-to-'
                                          'digital converter'))

    parent.append(odml.Section(name='Filter',
                               type='setup/daq',
                               definition='Information on the filter '
                                          'settings'))

    parent.append(odml.Section(name='Amplifier',
                               type='setup/daq',
                               definition='Information on the amplifier'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='Front-End Amplifier',
                                                 dtype=odml.DType.string),
                                definition='Type of amplifier'))

    parent.append(odml.Property(name='Model',
                                value=default_value(odml.DType.string),
                                definition='Model of amplifier'))

    parent.append(odml.Property(name='ChannelCount',
                                value=odml.Value(data=128,
                                                 dtype=odml.DType.int),
                                definition='Total number of analog input '
                                           'channels'))

    parent.append(odml.Property(name='UsedChannelCount',
                                value=odml.Value(data=96,
                                                 dtype=odml.DType.int),
                                definition='Number of analog input channels '
                                           'used'))

    parent.append(odml.Property(name='SavedIn',
                                value=odml.Value(data='ns6',
                                                 dtype=odml.DType.string),
                                definition='File in which data would be '
                                           'saved (if it is saved)'))

    parent.append(odml.Property(name='IsSaved',
                                value=default_value(odml.DType.boolean),
                                definition='States if signal is saved without '
                                           'further processing'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalAmplifier']['ADConverter']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='AIRange',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.float,
                                                  unit='mV')
                                       for data in [-8.192, 8.192]],
                                definition='Voltage range of analog input'))

    parent.append(odml.Property(name='MaxAIImpedance',
                                value=default_value(odml.DType.float,
                                                    unit=u"G\u2126"),
                                definition='Maximum impedance of analog '
                                           'input'))

    parent.append(odml.Property(name='Resolution',
                                value=odml.Value(data=0.25,
                                                 dtype=odml.DType.float,
                                                 unit=u'\u00B5V/bit'),
                                definition='Resolution of AD converter'))

    parent.append(odml.Property(name='DataBits',
                                value=odml.Value(data=16,
                                                 dtype=odml.DType.int,
                                                 unit='bit'),
                                definition='Width of digital data in bits'))

    parent.append(odml.Property(name='SamplingRate',
                                value=default_value(odml.DType.int,
                                                    unit='samples/sec'),
                                definition='Sampling rate of data'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalAmplifier']['Filter']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Type',
                                value=odml.Value(data='butterworth',
                                                 dtype=odml.DType.string),
                                definition='Filter type used'))

    parent.append(odml.Property(name='Causal',
                                value=odml.Value(data=True,
                                                 dtype=odml.DType.boolean),
                                definition='States if the used filter is '
                                           'causal'))

    parent.append(odml.Property(name='HighPassFreq',
                                value=odml.Value(data=0.3,
                                                 dtype=odml.DType.float,
                                                 unit='Hz'),
                                definition='High pass frequency used'))

    parent.append(odml.Property(name='LowPassFreq',
                                value=odml.Value(data=7.5,
                                                 dtype=odml.DType.float,
                                                 unit='kHz'),
                                definition='Low pass frequency used'))

    parent.append(odml.Property(name='HighPassOrder',
                                value=odml.Value(data=1,
                                                 dtype=odml.DType.int),
                                definition='High pass order used'))

    parent.append(odml.Property(name='LowPassOrder',
                                value=odml.Value(data=3,
                                                 dtype=odml.DType.int),
                                definition='Low pass order used'))

    # SET PARENT NODE
    parent = doc['Cerebus']['NeuralSignalAmplifier']['Amplifier']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Gain',
                                value=odml.Value(data=5000,
                                                 dtype=odml.DType.int),
                                definition='Gain used'))

    # SET PARENT NODE
    parent = doc['Cerebus']['ControlComputer']

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
    parent = doc['Cerebus']['ControlComputer']['ControlSoftware']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Name',
                                value=odml.Value(data='CentralSuite',
                                                 dtype=odml.DType.string),
                                definition='Name of software'))

    parent.append(odml.Property(name='Manufacturer',
                                value=odml.Value(data='Blackrock Microsystems',
                                                 dtype=odml.DType.string),
                                definition='Manufacturer of software'))

    parent.append(odml.Property(name='Version',
                                value=default_value(odml.DType.string),
                                definition='Software version'))

    # SET PARENT NODE
    parent = doc['UtahArray']

    # APPEND SUBSECTIONS
    parent.append(get_connector_subsection(connector_style=connector_style))

    parent.append(odml.Section(name='Array',
                               type='setup/daq/hardware',
                               definition='Information on the array used'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Owner',
                                value=odml.Value(data='Riken, Japan',
                                                 dtype=odml.DType.string),
                                definition='Owner of the device'))

    parent.append(odml.Property(name='SerialNo',
                                value=default_value(odml.DType.string),
                                definition='Serial number of the device'))

    parent.append(odml.Property(name='Manufacturer',
                                value=odml.Value(data='Blackrock Microsystems',
                                                 dtype=odml.DType.string),
                                definition='Manufacturer of the device'))

    parent.append(odml.Property(name='WiresMaterial',
                                value=odml.Value(data='Pt/Au',
                                                 dtype=odml.DType.string),
                                definition='Material of wires connecting '
                                           'array and connector'))

    parent.append(odml.Property(name='WiresDiameter',
                                value=odml.Value(data=25.0,
                                                 dtype=odml.DType.float,
                                                 unit=u'\u00B5m'),
                                definition='Diameter of each wire'))

    parent.append(odml.Property(name='WiresLength',
                                value=default_value(odml.DType.float,
                                                    unit='mm'),
                                definition='Length of each wire'))

    # SET PARENT NODE
    parent = doc['UtahArray']['Array']

    # APPEND SUBSECTIONS
    grid_section_names = get_grid_section_names(grid_no)
    for sname in grid_section_names:
        parent.append(odml.Section(name=sname,
                                   type='setup/daq/hardware',
                                   definition='Information on array grid'))

    elec_section_names = get_elec_section_names(electrode_ids)
    for sname in elec_section_names:
        parent.append(odml.Section(name=sname,
                                   type='setup/daq/hardware',
                                   definition='Information on electrode'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='GridCount',
                                value=odml.Value(data=grid_no,
                                                 dtype=odml.DType.int),
                                definition='Number of array grids'))

    parent.append(odml.Property(name='ReferenceCount',
                                value=odml.Value(data=2,
                                                 dtype=odml.DType.int),
                                definition='Number of reference wires'))

    parent.append(odml.Property(name='GroundCount',
                                value=odml.Value(data=1,
                                                 dtype=odml.DType.int),
                                definition='Number of ground wires'))

    parent.append(odml.Property(name='ElectrodeCount',
                                value=odml.Value(data=100,
                                                 dtype=odml.DType.int),
                                definition='Total number of electrodes'))

    parent.append(odml.Property(name='ActiveElectrodeCount',
                                value=odml.Value(data=96,
                                                 dtype=odml.DType.int),
                                definition='Number of connected/active '
                                           'electrodes'))

    parent.append(odml.Property(name='UsedElectrodeCount',
                                value=odml.Value(data=len(electrode_ids),
                                                 dtype=odml.DType.int),
                                definition='Number of used/saved electrodes'))

    for sname in grid_section_names:
        # SET PARENT NODE
        parent = doc['UtahArray']['Array'][sname]

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='ID',
                                    value=default_value(odml.DType.int)
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data=int(sname[-2:]),
                                                    dtype=odml.DType.int),
                                    definition='Grid ID'))

        parent.append(odml.Property(name='ElectrodeGeometry',
                                    value=default_value(odml.DType.string)
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data='flat',
                                                    dtype=odml.DType.string),
                                    definition='Geometry of the electrode '
                                               'lengths'))

        parent.append(odml.Property(name='ElectrodeMetal',
                                    value=default_value(odml.DType.string)
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data='iridium oxide',
                                                    dtype=odml.DType.string),
                                    definition='Electrode site metal option'))

        parent.append(odml.Property(name='ElectrodePitch',
                                    value=default_value(odml.DType.float,
                                                        unit=u'\u00B5m')
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data=400.0,
                                                    dtype=odml.DType.float,
                                                    unit=u'\u00B5m'),
                                    definition='Inter-electrode-distance'))

        parent.append(odml.Property(name='Insulation',
                                    value=odml.Value(data='Parylene-C',
                                                     dtype=odml.DType.string),
                                    definition='Insulation material'))

        parent.append(odml.Property(name='GridRows',
                                    value=default_value(odml.DType.int)
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data=10,
                                                    dtype=odml.DType.int),
                                    definition='Number of electrode rows.'))

        parent.append(odml.Property(name='GridColumns',
                                    value=default_value(odml.DType.int)
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data=10,
                                                    dtype=odml.DType.int),
                                    definition='Number of electrode columns'))

        parent.append(odml.Property(name='GridWidth',
                                    value=default_value(odml.DType.float,
                                                        unit='mm')
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data=4.0,
                                                    dtype=odml.DType.float,
                                                    unit='mm'),
                                    definition='Width of array grid'))

        parent.append(odml.Property(name='GridLength',
                                    value=default_value(odml.DType.float,
                                                        unit='mm')
                                    if sname[-2:] == 'XX'
                                    else odml.Value(data=4.0,
                                                    dtype=odml.DType.float,
                                                    unit='mm'),
                                    definition='Length of array grid'))

        parent.append(odml.Property(name='ElectrodeCount',
                                    value=odml.Value(data=100,
                                                     dtype=odml.DType.int),
                                    definition='Total number of electrodes'))

        parent.append(odml.Property(name='ActiveElectrodeCount',
                                    value=odml.Value(data=len(electrode_ids),
                                                     dtype=odml.DType.int),
                                    definition='Number of connected/active '
                                               'electrodes'))

    for sname in elec_section_names:
        # SET PARENT NODE
        parent = doc['UtahArray']['Array'][sname]

        # APPEND SUBSECTIONS
        if spike_sorted:
            parent.append(odml.Section(name='OfflineSpikeSorting',
                                       type='setup/daq/hardware',
                                       definition='Information on the spike '
                                                  'sorting performed offline, '
                                                  'after the experiment'))

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='ID',
                                    value=default_value(odml.DType.int)
                                    if sname[-3:] == 'XXX'
                                    else odml.Value(data=int(sname[-3:]),
                                                    dtype=odml.DType.int),
                                    definition='Electrode ID'))

        parent.append(odml.Property(name='GridID',
                                    value=default_value(odml.DType.int),
                                    definition='Grid ID'))

        parent.append(odml.Property(name='BankID',
                                    value=default_value(odml.DType.int),
                                    definition='Bank ID'))

        parent.append(odml.Property(name='PinID',
                                    value=default_value(odml.DType.int),
                                    definition='Pin ID'))

        parent.append(odml.Property(name='IDca',
                                    value=default_value(odml.DType.int),
                                    definition='Connector aligned ID (array '
                                    'head view, wiring to connector on right '
                                    'side, count order from left to right and '
                                    ' bottom to top)'))

        parent.append(odml.Property(name='IDba',
                                    value=default_value(odml.DType.int),
                                    definition='Brain aligned ID (central '
                                    'sulcus at top)'))

        parent.append(odml.Property(name='Impedance',
                                    value=default_value(odml.DType.float,
                                                        unit=u'k\u2126'),
                                    definition='Pre-implantation impedance'))

        parent.append(odml.Property(name='AutoImpedance',
                                    value=default_value(odml.DType.float,
                                                        unit=u'k\u2126'),
                                    definition='Auto impedance during '
                                               'implantation'))

        parent.append(odml.Property(name='Length',
                                    value=default_value(odml.DType.float,
                                                        unit='mm')
                                    if sname[-3:] == 'XXX'
                                    else odml.Value(data=1.5,
                                                    dtype=odml.DType.float,
                                                    unit='mm'),
                                    definition='Length'))

        if spike_sorted:
            # SET PARENT NODE
            parent = doc['UtahArray']['Array'][sname]['OfflineSpikeSorting']

            # APPEND SUBSECTIONS

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='SUAIDs',
                                        value=default_value(odml.DType.int),
                                        definition='IDs of single units'))

            parent.append(odml.Property(name='MUAID',
                                        value=odml.Value(data=0,
                                                         dtype=odml.DType.int),
                                        definition='IDs of multi units'))

            parent.append(odml.Property(name='NoiseIDs',
                                        value=odml.Value(data=255,
                                                         dtype=odml.DType.int),
                                        definition='IDs for rejected units'))

            parent.append(odml.Property(name='SpikeThreshold',
                                        value=default_value(odml.DType.float,
                                                            unit='mV'),
                                        definition='Threshold used for spike '
                                                   'detection'))

    # SET PARENT NODE
    parent = doc['Headstage']

    # APPEND SUBSECTIONS

    headstage_sets = get_headstage_sets(used_headstage)
    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Model',
                                value=headstage_sets['Model'],
                                definition='Model of device'))

    parent.append(odml.Property(name='Type',
                                value=headstage_sets['Type'],
                                definition='Type of device (digital or '
                                           'analog)'))

    parent.append(odml.Property(name='Gain',
                                value=odml.Value(data=1,
                                                 dtype=odml.DType.int)
                                if headstage_sets
                                else default_value(odml.DType.int),
                                definition='Gain used'))

    parent.append(odml.Property(name='InFrom',
                                value=headstage_sets['InFrom'],
                                definition='Devices from which inputs are '
                                           'received'))

    parent.append(odml.Property(name='OutTo',
                                value=headstage_sets['OutTo'],
                                definition='Target device of each output'))

    parent.append(odml.Property(name='StimSwitch',
                                value=default_value(odml.DType.boolean),
                                definition='States if a stimulation switch '
                                           'version of the device was used'))

    return doc
