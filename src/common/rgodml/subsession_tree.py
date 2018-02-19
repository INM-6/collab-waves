# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 10:09:15 2014

@author: zehl
"""
import odml
from rgodml.odml_io import default_value


def get_trial_section_names(trialids=[]):
    """
    Returns a list of trial Section names

    Args:
        trialids (list of int, optional):
            List of trial ids to use for trial Section names

    Returns:
        (list of strings):
            List of trial Section names

    Note:
        Trial Section names are formatted as 'Trial_ID' where ID is given as a
        three digit id created from the given trial id list (e.g. ['Trial_001',
        'Trial_003'] for trialids = [1, 3])

        If len(trialids) = 0 the list ['Trial_XXX'] is returned
    """
    if len(trialids) == 0:
        section_names = ['Trial_XXX']

    elif len(trialids) > 0:
        section_names = ['Trial_%03i' % i for i in trialids]

    else:
        raise ValueError('Trial IDs must be a list of integers')

    return section_names


def get_task_subsection(tasktype, trialids=[]):
    """
    Constructs a Section describing the given task type including Properties,
    Values and (possibly default) data

    Args:
        tasktype (string):
            Type of performed task (OneCue, TwoCues, Observation, Sleep, or
            Mapping)
        trialids (list of int):
            trial ids (only possible for OneCue, TwoCues, Observation)

    Returns:
        (odml.section.BaseSection)
            Section describing the specified task type

    Note:
        Valid task types: OneCue, TwoCues, Observation, Sleep, Mapping
    """
    # CREATE MAIN SECTION
    main = odml.Section(name='TaskSettings',
                        type='subsession/task',
                        definition='Information on task settings')

    if tasktype in ['OneCue', 'TwoCues', 'Observation']:
        # APPEND SUBSECTIONS
        main.append(odml.Section(name='TrialTypeSettings',
                                 type='subsession/task/trial',
                                 definition='Information on trial type '
                                            'parameters'))

        trial_section_names = get_trial_section_names(trialids)
        for sname in trial_section_names:
            main.append(odml.Section(name=sname,
                                     type='subsession/task/trial',
                                     definition='Information on trial '
                                                'parameters'))

        # APPEND PROPERTIES WITH DEFAULT VALUES
        main.append(odml.Property(name='TotalTrialCount',
                                  value=default_value(odml.DType.int),
                                  definition='Total number of trials'))

        main.append(odml.Property(name='CorrectTrialCount',
                                  value=default_value(odml.DType.int),
                                  definition='Number of correctly performed '
                                             'trials'))

        main.append(odml.Property(name='GripErrorTrialCount',
                                  value=default_value(odml.DType.int),
                                  definition='Number of error trials'))

        main.append(odml.Property(name='StandardSettings',
                                  value=default_value(odml.DType.boolean),
                                  definition='States if subsession was '
                                             'performed under standard '
                                             'settings'))
        # SET PARENT NODE
        parent = main['TrialTypeSettings']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Condition',
                                    value=default_value(odml.DType.int),
                                    definition='Trial type combination of '
                                               'subsession was condition'))

        parent.append(odml.Property(name='OrderForce',
                                    value=default_value(odml.DType.string),
                                    definition='Trial order of force '
                                               'conditions (random or block)'))

        parent.append(odml.Property(name='OrderGrip',
                                    value=default_value(odml.DType.string),
                                    definition='Trial order of grip '
                                               'conditions (random or block)'))

        parent.append(odml.Property(name='BlockSize',
                                    value=default_value(odml.DType.int),
                                    definition='Order of trial type '
                                               'combination was of block '
                                               'size'))

        # SET PARENT NODE
        for sname in trial_section_names:
            parent = main[sname]

            # APPEND SUBSECTIONS
            parent.append(odml.Section(name='DigitalEvents',
                                       type='subsession/task/trial/events',
                                       definition='Information on digital '
                                                  'trial events'))

            parent.append(odml.Section(name='AnalogEvents',
                                       type='subsession/task/trial/events',
                                       definition='Information on analogue '
                                                  'trial events'))

            parent.append(odml.Section(name='Periods',
                                       type='subsession/task/trial/periods',
                                       definition='Information on trial '
                                                  'periods'))

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='ID',
                                        value=default_value(odml.DType.int),
                                        definition='Trial ID'))

            parent.append(odml.Property(name='PerformanceCode',
                                        value=default_value(odml.DType.int),
                                        definition='States if this was an '
                                                   'error or a correct trial'))

            parent.append(odml.Property(name='TrialType',
                                        value=default_value(odml.DType.int),
                                        definition='Trial type of trial'))

            parent.append(odml.Property(name='ForceType',
                                        value=default_value(odml.DType.string),
                                        definition='Force type measured in '
                                                   'trial'))

            # SET PARENT NODE
            parent = main[sname]['DigitalEvents']

            # APPEND SUBSECTIONS

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='TS',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when trial start is '
                                                   'activated'))

            parent.append(odml.Property(name='FP-ON',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when fixation point '
                                                   'is activated'))

            parent.append(odml.Property(name='CUE-ON',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when (first) CUE is '
                                                   'activated '))

            parent.append(odml.Property(name='CUE-OFF',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when (first) CUE is '
                                                   'deactivated '))

            parent.append(odml.Property(name='GO-ON',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when (second cue) GO '
                                                   'is activated '))

            parent.append(odml.Property(name='SR',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when switch is '
                                                   'released'))

            parent.append(odml.Property(name='GO-OFF',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when (second cue) GO '
                                                   'is deactivated '))

            parent.append(odml.Property(name='RW',
                                        value=default_value(odml.DType.int,
                                                            unit='1./30000*s'),
                                        definition='Time when reward is '
                                                   'given'))

            # SET PARENT NODE
            parent = main[sname]['AnalogEvents']

            # APPEND SUBSECTIONS
            parent.append(odml.Section(name='GripForceSignals',
                                       type='subsession/task/trial/events',
                                       definition='Information on analogue '
                                                  'events detected with grip '
                                                  'force signals'))

            parent.append(odml.Section(name='DisplacementSignal',
                                       type='subsession/task/trial/events',
                                       definition='Information on analogue '
                                                  'events detected with '
                                                  'object displacement '
                                                  'signal'))

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='UsedForceSensor',
                                        value=default_value(odml.DType.int),
                                        definition='Channel ID index of force '
                                                   'sensor used in grip force '
                                                   'detection'))

            parent.append(odml.Property(name='ManuallyControlled',
                                        value=default_value(
                                            odml.DType.boolean),
                                        definition='States if analogue event '
                                                   'detection of this trial '
                                                   'was manually controlled'))

            # SET PARENT NODE
            parent = main[sname]['AnalogEvents']['GripForceSignals']

            # APPEND SUBSECTIONS

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='OT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when target object '
                                                   'is touched'))

            parent.append(odml.Property(name='HS',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when holding '
                                                   'position of target object '
                                                   'is reached and holding '
                                                   'period starts'))

            parent.append(odml.Property(name='OR',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when target object '
                                                   'is released'))

            parent.append(odml.Property(name='BTB',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when object moved '
                                                   'back to its baseline'))

            # SET PARENT NODE
            parent = main[sname]['AnalogEvents']['DisplacementSignal']

            # APPEND SUBSECTIONS

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='OT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when target object '
                                                   'is touched'))

            parent.append(odml.Property(name='HS',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when holding '
                                                   'position of target object '
                                                   'is reached and holding '
                                                   'period starts'))

            parent.append(odml.Property(name='OR',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when target object '
                                                   'is released'))

            parent.append(odml.Property(name='BTB',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time when object moved '
                                                   'back to its baseline'))

            # SET PARENT NODE
            parent = main[sname]['Periods']

            # APPEND SUBSECTIONS
            parent.append(odml.Section(name='GripForceSignals',
                                       type='subsession/task/trial/events',
                                       definition='Information on periods '
                                                  'calculated with analogue '
                                                  'events detected with grip '
                                                  'force signals'))

            parent.append(odml.Section(name='DisplacementSignal',
                                       type='subsession/task/trial/events',
                                       definition='Information on periods '
                                                  'calculated with analogue '
                                                  'events detected with '
                                                  'object displacement '
                                                  'signal'))

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='RT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Reaction time (period '
                                                   'between GO and SR)'))

            parent.append(odml.Property(name='ITI',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Inter-trial-interval '
                                                   '(period until next '
                                                   'trial)'))

            # SET PARENT NODE
            parent = main[sname]['Periods']['GripForceSignals']

            # APPEND SUBSECTIONS

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='RGT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Reach-to-grasp time '
                                                   '(period between SR and '
                                                   'OT)'))

            parent.append(odml.Property(name='PT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Pull time (period between '
                                                   'OT and HS)'))

            parent.append(odml.Property(name='MT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Movement time (period '
                                                   'between GO-ON and HS)'))

            parent.append(odml.Property(name='HT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Hold time (period between '
                                                   'HS and RW)'))

            parent.append(odml.Property(name='ORT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time until object is '
                                                   'released (period between '
                                                   'RW an OR '))

            # SET PARENT NODE
            parent = main[sname]['Periods']['DisplacementSignal']

            # APPEND SUBSECTIONS

            # APPEND PROPERTIES WITH DEFAULT VALUES
            parent.append(odml.Property(name='RGT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Reach-to-grasp time '
                                                   '(period between SR and '
                                                   'OT)'))

            parent.append(odml.Property(name='PT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Pull time (period between '
                                                   'OT and HS)'))

            parent.append(odml.Property(name='MT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Movement time (period '
                                                   'between GO-ON and HS)'))

            parent.append(odml.Property(name='HT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Hold time (period between '
                                                   'HS and RW)'))

            parent.append(odml.Property(name='ORT',
                                        value=default_value(odml.DType.int,
                                                            unit='ms'),
                                        definition='Time until object is '
                                                   'released (period between '
                                                   'RW an OR '))

    elif tasktype == 'Sleep':
        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        main.append(odml.Property(name='SleepStart',
                                  value=default_value(odml.DType.int,
                                                      unit='samples'),
                                  definition='Start time(s) of sleep '
                                             'interval(s)'))

        main.append(odml.Property(name='SleepEnd',
                                  value=default_value(odml.DType.int,
                                                      unit='samples'),
                                  definition='End time(s) of sleep '
                                             'interval(s)'))

    elif tasktype == 'Mapping':
        # APPEND SUBSECTIONS
        main.append(odml.Section(name='Shoulder',
                                 type='subsession/task/RF',
                                 definition='Information on shoulder mapping'))

        main.append(odml.Section(name='Elbow',
                                 type='subsession/task/RF',
                                 definition='Information on elbow mapping'))

        main.append(odml.Section(name='Wrist',
                                 type='subsession/task/RF',
                                 definition='Information on wrist mapping'))

        main.append(odml.Section(name='Thumb',
                                 type='subsession/task/RF',
                                 definition='Information on wrist mapping'))

        main.append(odml.Section(name='Index',
                                 type='subsession/task/RF',
                                 definition='Information on wrist mapping'))

        main.append(odml.Section(name='Digit2to5',
                                 type='subsession/task/RF',
                                 definition='Information on wrist mapping'))

        main.append(odml.Section(name='Digit5',
                                 type='subsession/task/RF',
                                 definition='Information on wrist mapping'))

        # APPEND PROPERTIES WITH DEFAULT VALUES
        main.append(odml.Property(name='PassivePG',
                                  value=default_value(odml.DType.int,
                                                      unit='samples'),
                                  definition='Time interval where a '
                                             'passively precision grip was '
                                             'mimicked'))

        main.append(odml.Property(name='UlnarDeviation',
                                  value=default_value(odml.DType.int,
                                                      unit='samples'),
                                  definition='Time interval where Ulnar was '
                                             'deviated'))

        main.append(odml.Property(name='Digit2to3Extension',
                                  value=default_value(odml.DType.int,
                                                      unit='samples'),
                                  definition='Time interval where digits 2 '
                                             'and 3 were extended'))

        main.append(odml.Property(name='TouchPalm',
                                  value=default_value(odml.DType.int,
                                                      unit='samples'),
                                  definition='Time interval where palm was '
                                             'touched'))

        main.append(odml.Property(name='TouchCheek',
                                  value=default_value(odml.DType.int,
                                                      unit='samples'),
                                  definition='Time interval where cheek was '
                                             'touched'))

        # SET PARENT NODE
        parent = main['Shoulder']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='MoveForward',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was moved forward'))

        parent.append(odml.Property(name='Lowering',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was lowered down'))

        # SET PARENT NODE
        parent = main['Elbow']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Flexion',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was flexed'))

        parent.append(odml.Property(name='Extension',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was extended'))

        # SET PARENT NODE
        parent = main['Wrist']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Flexion',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was flexed'))

        parent.append(odml.Property(name='Pronation',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was pronated'))

        # SET PARENT NODE
        parent = main['Thumb']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Abduction',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was abducted'))

        parent.append(odml.Property(name='Adduction',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was adducted'))

        parent.append(odml.Property(name='Extension',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which body '
                                               'part was extended'))

        parent.append(odml.Property(name='TouchEnd',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which end of '
                                               'body part was touched'))

        parent.append(odml.Property(name='TouchInterior',
                                    value=default_value(odml.DType.int,
                                                        unit='samples'),
                                    definition='Time interval in which '
                                               'interior body part was '
                                               'touched'))

        # SET PARENT NODE
        parent = main['Index']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Abduction',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which body part was '
                      'abducted'))

        parent.append(odml.Property(name='Flexion',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which body part was '
                      'flexed'))

        parent.append(odml.Property(name='Extension',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which body part was '
                      'extended'))

        parent.append(odml.Property(name='TouchInteriorFace',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which interior body part '
                      'was touched'))

        # SET PARENT NODE
        parent = main['Digit2to5']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Flexion',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which body part was '
                      'flexed'))

        parent.append(odml.Property(name='TouchInteriorFace',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which interior body part '
                      'was touched'))

        # SET PARENT NODE
        parent = main['Digit5']

        # APPEND SUBSECTIONS

        # APPEND PROPERTIES WITH DEFAULT VALUES
        parent.append(odml.Property(name='Flexion',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which body part was '
                      'flexed'))

        parent.append(odml.Property(name='TouchLateralFace',
                      value=default_value(odml.DType.int, unit='samples'),
                      definition='Time interval in which lateral body part '
                      'was touched'))

    return main


def rg_subsession_tree(tasktype=None, trialids=[]):
    """
    Builds an odML Document for reach-to-grasp subsession metadata

    Returns:
        (odml.doc.BaseDocument)
            odML Document with default reach-to-grasp subsession metadata
    """

    # CREATE A DOCUMENT
    doc = odml.Document(version='0.1')

    # APPEND MAIN SECTION
    doc.append(odml.Section(name='Subsession',
                            type='subsession',
                            definition='Information on the subsession'))

    # SET PARENT NODE
    parent = doc['Subsession']

    # APPEND SUBSECTIONS
    if tasktype:
        parent.append(get_task_subsection(tasktype=tasktype,
                                          trialids=trialids))

    # APPEND PROPERTIES WITH DEFAULT VALUES
    parent.append(odml.Property(name='Session',
                                value=default_value(odml.DType.string),
                                definition='Session name'))

    parent.append(odml.Property(name='Subsession',
                                value=default_value(odml.DType.string),
                                definition='Subsession name'))

    parent.append(odml.Property(name='Date',
                                value=default_value(odml.DType.date),
                                definition='Date of session'))

    parent.append(odml.Property(name='Weekday',
                                value=default_value(odml.DType.string),
                                definition='Weekday of session'))

    parent.append(odml.Property(name='Time',
                                value=default_value(odml.DType.time),
                                definition='Time of subsession'))

    parent.append(odml.Property(name='Duration',
                                value=default_value(odml.DType.int,
                                                    unit='sec'),
                                definition='Duration of subsession'))

    parent.append(odml.Property(name='IsSpikeSorted',
                                value=default_value(odml.DType.boolean),
                                definition='States if a spike sorting was '
                                           'performed for this subsession '))

    value = default_value(odml.DType.string)
    if tasktype:
        value.data = tasktype
    parent.append(odml.Property(name='TaskType',
                                value=value,
                                definition='Performed task type'))

    parent.append(odml.Property(name='Noisy',
                                value=default_value(odml.DType.boolean),
                                definition='States if subsession is noisy'))

    parent.append(odml.Property(name='Comment',
                                value=default_value(odml.DType.string),
                                definition='Comment about the subsession'))

    return doc
