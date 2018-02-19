# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 12:23:31 2014

@author: zehl
"""
import odml
from rgodml.odml_io import default_value


def rej_tree(rej_names=[]):
    """
    Builds an odML Document for reach-to-grasp preprocessing metadata

    Returns:
        (odml.doc.BaseDocument)
            odML Document with default reach-to-grasp preprocessing metadata
    """

    # CREATE A DOCUMENT
    doc = odml.Document(version='0.1')

    # APPEND MAIN SECTION
    doc.append(odml.Section(name='RejectionsLFP',
                            type='preprocessing',
                            definition='Information on rejection of '
                                       'experimental LFP signals'))

    # SET PARENT NODE
    parent = doc['RejectionsLFP']

    # APPEND SUBSECTIONS
    for rejn in rej_names:
        parent.append(odml.Section(name=rejn,
                                   type='subject/preparation',
                                   definition='Information on rejection '
                                              'performed on this LFP band'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='LFPBands',
                                value=[odml.Value(data=data,
                                                  dtype=odml.DType.string)
                                       for data in rej_names]
                                if len(rej_names) > 0
                                else default_value(odml.DType.string),
                                definition='Information on which LFP bands '
                                           'the rejection was performed'))

    if len(rej_names) > 0:
        for rejn in rej_names:
            # SET PARENT NODE
            parent = doc['RejectionsLFP'][rejn]

            # APPEND SUBSECTIONS

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='Highcut',
                                        value=default_value(odml.DType.float),
                                        definition='High cut frequency of '
                                                   'rejection filter'))

            parent.append(odml.Property(name='Lowcut',
                                        value=default_value(odml.DType.float),
                                        definition='Low cut frequency of '
                                                   'rejection filter'))

            parent.append(odml.Property(name='Order',
                                        value=default_value(odml.DType.int),
                                        definition='Order of rejection filter'))

            parent.append(odml.Property(name='RejElectrodes',
                                        value=default_value(odml.DType.int),
                                        definition='ID list of rejected '
                                                   'electrodes'))

            parent.append(odml.Property(name='RejTrials',
                                        value=default_value(odml.DType.int),
                                        definition='ID list of rejected '
                                                   'trials'))

    return doc