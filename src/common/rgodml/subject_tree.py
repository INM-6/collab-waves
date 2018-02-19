# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 17:48:31 2014

@author: zehl
"""
import odml
from rgodml.odml_io import default_value


def rg_subject_tree():
    """
    Builds an odML Document for reach-to-grasp subject metadata

    Returns:
        (odml.doc.BaseDocument)
            odML Document with default reach-to-grasp subject metadata
    """

    # CREATE A DOCUMENT
    doc = odml.Document(version='0.1')

    # APPEND MAIN SECTION
    doc.append(odml.Section(name='Subject',
                            type='subject',
                            definition='Information on the investigated '
                                       'experimental subject (animal or '
                                       'person)'))

    # SET PARENT NODE
    parent = doc['Subject']

    # APPEND SUBSECTIONS
    parent.append(odml.Section(name='Training',
                               type='subject/preparation',
                               definition='Information on the training given '
                                          'to subject'))

    parent.append(odml.Section(name='ArrayImplant',
                               type='subject/preparation',
                               definition='Information on the array implant '
                                          'performed on subject'))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Species',
                                value=default_value(odml.DType.string),
                                definition='Binomial species name (genus, '
                                           'species within genus)'))

    parent.append(odml.Property(name='TrivialName',
                                value=default_value(odml.DType.string),
                                definition='Commonly used (trivial) species '
                                           'name'))

    parent.append(odml.Property(name='GivenName',
                                value=default_value(odml.DType.string),
                                definition='Given name'))

    parent.append(odml.Property(name='Identifier',
                                value=default_value(odml.DType.string),
                                definition='Identifier used in file names'))

    parent.append(odml.Property(name='Gender',
                                value=default_value(odml.DType.string),
                                definition='Gender (male or female)'))

    parent.append(odml.Property(name='Birthday',
                                value=default_value(odml.DType.date),
                                definition='Date of birth (yyyy-mm-dd)'))

    parent.append(odml.Property(name='ActiveHand',
                                value=default_value(odml.DType.string),
                                definition='Trained hand (left and/or right)'))

    parent.append(odml.Property(name='Disabilities',
                                value=default_value(odml.DType.text),
                                definition='Comment on existing disabilities'))

    parent.append(odml.Property(name='Character',
                                value=default_value(odml.DType.text),
                                definition='Comment on the general character/'
                                           'personality'))

    # SET PARENT NODE
    parent = doc['Subject']['Training']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Coach',
                                value=default_value(odml.DType.string),
                                definition='Full name(s) of person(s) '
                                           'responsible for the training'))

    parent.append(odml.Property(name='Start',
                                value=default_value(odml.DType.date),
                                definition='Start date of training'))

    parent.append(odml.Property(name='End',
                                value=default_value(odml.DType.date),
                                definition='End date of training'))

    parent.append(odml.Property(name='Protocol',
                                value=default_value(odml.DType.url),
                                definition='Training protocol URL'))

    parent.append(odml.Property(name='Comment',
                                value=default_value(odml.DType.string),
                                definition='Comment on training'))

    # SET PARENT NODE
    parent = doc['Subject']['ArrayImplant']

    # APPEND SUBSECTIONS

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Date',
                                value=default_value(odml.DType.date),
                                definition='Date of the surgery'))

    parent.append(odml.Property(name='Surgeon',
                                value=default_value(odml.DType.string),
                                definition='Full name(s) of the person(s) '
                                           'responsible for the surgery'))

    parent.append(odml.Property(name='Protocol',
                                value=default_value(odml.DType.url),
                                definition='Surgery protocol URL'))

    parent.append(odml.Property(name='Pictures',
                                value=default_value(odml.DType.url),
                                definition='URL of pictures taken during the '
                                           'surgery'))

    parent.append(odml.Property(name='Hemisphere',
                                value=default_value(odml.DType.string),
                                definition='Hemisphere (left and/or right) in '
                                           'which surgery was performed'))

    parent.append(odml.Property(name='Comment',
                                value=default_value(odml.DType.string),
                                definition='Comment on the surgery'))

    return doc
