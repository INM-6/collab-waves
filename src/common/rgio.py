# coding=utf-8
'''
Reach-grasp IO functions

This module can be used to load data from the reach grasp experiments via the
Blackrock loading IO of neo.

Authors: Michael Denker, Lyuba Zehl
'''

import copy
import glob
import itertools
import os
import sys
import warnings

import numpy as np
import scipy.signal
import odml.tools
import quantities as pq

import neo
from neo.io.blackrockio import BlackrockIO
import rgodml.metadata_io


class ReachGraspIO(BlackrockIO):
    """
    Derived class to handle a managing a Blackrock session recording from the
    reach to grasp experiments.

    Args:
        filename (string):
            Name of a reach-grasp Blackrock file to associate with (i.e.,
            the subsession). The file extension should not be included.
        sorting_version (string):
            Postfix to the subsession name to identify the nev file to use
            if several nev files exist. (example: filename='a' and
            sorting_version='-test-05' will load 'a-test-05.nev'). If set to
            None, the best version will be used. The best version is one
            where a .txt file exists, the filename includes "test", and the
            last number in the nev file name is highest. To use the original
            Blackrock nev, set sorting_version="". Default: None.
        metadata_dir (string):
            Directory where the metadata information (including the .odml/.json
            files for the subsession) is stored. This directory is typically
            called "MetaDataMonkeyname" and must contain the 'odMLbuild' and
            'source' subdirectories. The filename of the odML file must equal
            the filename of the subsession and it must reside in the
            'odMLbuild' subdirectory. If no odML file is found in the odMLbuild
            directory, the loader will attempt to extract some basic metadata
            from the 'source' subdirectory of metadata_dir. If None is
            specified, no odML or metadata information is loaded.
            Default: None
        print_diagnostic (boolean):
            If true, the routine will output diagnostic information. This
            option will likely be removed in the future!
            Default: False

    Attributes:
        txt_fileprefix (str):
            File name of the chosen .txt file (without extension).
        num_total_trials (int):
            Number of trials in the file, where a trial is defined by a reward
            event.
        num_trials (int):
            Number of correct trials, i.e., those trials which have a
            correct order of the events from trial start to reward.
        is_one_cue_task (bool):
            True if this subsession contains only one cue. In this case the
            delay cue is set to the go cue.
        is_two_cues_task (bool):
            True if this subsession contains only two different cues. In this
            case the delay cue is of the opposite type then the go cue. The
            cues can be either of type force or of type grip.
        is_sorted (bool):
            True if a genuine sorting (i.e., a txt file) has been loaded to
            distinguish SUA and MUA.
        trial_data (numpy.array):
            A num_trialsx10 matrix, where each row corresponds to one trial.
            Columns:
                 0: Consecutive trial index determined by trial start triggers,
                         i.e., including numbered error trials (trial_id)
                 1: Depricated old trial index determined by reward triggers,
                         i.e., excluding most error trials (trial_id_depr)
                 2: Performance code of each trial determined by ability to
                         detect the trial events, i.e. 127 for a correct
                         trial with all events occurring, see performance_str
                         below for codes, (performance_code)
                 3: trial type, see trial_types below for codes (trial_type)
                 4: det_force, as detected from force signals (1=LF,2=HF, or
                         0=not available)
                 5-14: Time of events in samples (starting at 0):
                     3: TS (trial start, self initiated by monkey)
                     4: FP-ON (fixation cue on)
                     5: CUE-ON (first/delay cue on)
                     6: CUE-OFF (first/delay cue off)
                     7: GO-ON (second/go cue on)
                     8: GO-OFF (second/go cue off)
                     9: SR (switch release, start of reaching)
                    10: RW (reward)
                    11: OT (object touch, end of reaching,
                            start of pull and hold)
                    12: OR (object release, end of holding after reward)
        trial_data_idx (dict):
            Index into trial_data (value) by column name (key) (see trial_data
            above).
        trial_data_str (dict):
            Name of trial_data column (value) by index (key) (see trial_data
            above).
        trial_types (list):
            List of all occurring trial type codes in the subsession. Possible
            trial types:
                85:  SGSG (first cue: side grip, second cue: side grip)
                86:  SGHF (first cue: side grip, second cue: high force)
                89:  SGLF (first cue: side grip, second cue: low force)
                101: HFSG (first cue: high force, second cue: side grip)
                102: HFHF (first cue: high force, second cue: high force)
                106: HFPG (first cue: high force, second cue: precision grip)
                149: LFSG (first cue: low force, second cue: side grip)
                153: LFLF (first cue: low force, second cue: low force)
                154: LFPG (first cue: low force, second cue: precision grip)
                166: PGHF (first cue: precision grip, second cue: high force)
                169: PGLF (first cue: precision grip, second cue: low force)
                170: PGPG (first cue: precision grip, second cue: precision
                        grip)
        trial_type_str (dict):
            Short-hand notation of the trial type (value) for the integer trial
            type code (keys). (see trial_types)
        trial_type_codes (dict):
            Integer trial type code (value) for the short-hand notation of the
            trial type (key). (see trial_types)
        performance_code_str (dict):
            Short-hand notation of performance code (value) for the integer
            code of trial performance (keys):
                127: correct_trials  (all events occur correctly)
                64:  SR<FP-ON        (SR before FP-ON)
                96:  SR<CUE-ON       (SR before CUE-ON)
                112: SR<CUE-OFF      (SR before CUE-OFF)
                120: SR<GO           (SR before GO-ON event)
                124: SR_missing      (no SR is detected)
                126: grip_error      (all events occur correctly, but RW)
        performance_code_idx (dict):
            Integer code of trial performance (value) for the short-hand
            notation of performance code (key). (see performance_code_str)
        condition (int):
            Code indicating the condition. Possible conditions:
               0:
                 No trials, or condition not conclusive from file
               4 types (two_cues_task):
                 1: all grip-first trial types with two different cues
                 2: all force-first trial types with two different cues
               2 types (two_cues_task):
                 11: grip-first, but only LF types
                 12: grip-first, but only HF types
                 13: grip-first, but only SG types
                 14: grip-first, but only PG types
               2 types (two_cues_task):
                 21: force-first, but only LF types
                 22: force-first, but only HF types
                 23: force-first, but only SG types
                 24: force-first, but only PG types
               1 type (two_cues_task):
                 131: grip-first, but only SGLF type
                 132: grip-first, but only SGHF type
                 141: grip-first, but only PGLF type
                 142: grip-first, but only PGHF type
                 213: force-first, but only LFSG type
                 214: force-first, but only LFPG type
                 223: force-first, but only HFSG type
                 224: force-first, but only HFPG type
               1 type (one_cue_task):
                 133: SGSG, only grip info, force unknown
                 144: PGPG, only grip info, force unknown
                 211: LFLF, only force info, grip unknown
                 222: HFHF, only force info, grip unknown
        condition_trial_types_codes (dict):
            For each condition code (dict. keys), returns a list of trial type
            codes that occur in the condition.
            Example: condition_trial_types_codes[11] -> [89, 169]
        condition_trial_types_str (dict):
            For each condition code (dict. keys), returns a list of trial type
            shorthand strings that occur in the condition.
            Example: condition_trial_types_str[11] -> ['SGLF', 'PGLF']
        trial_events (list):
            List of all trial events as shorthand strings. (see trial_data)
        trial_events_str (dict):
            A dictionary of reach-grasp event shorthands as strings (values)
            and their indices in trial_data (keys). (see trial_data)
        trial_events_codes (dict):
            A dictionary of reach-grasp event indices in trial_data (values)
            and their shorthands as strings (keys). (see events_str,
            trial_data)
        force_channels (list):
            Indices of non-LFP analog channels
               137: analog force signal 1 for OT.OR calculation for PG trials
               138: analog force signal 1 for OT/OR calculation for SG trials
               139: analog force signal 2 for OT/OR calculation for PG trials
               140: analog force signal 2 for OT/OR calculation for SG trials
        force_channels_idx (dict):
            Dictionary of non-LFP analog channels shorthands (keys) with their
            corresponding channel indices (values). (see force_channels)
        force_channels_str (dict):
            Dictionary of non-LFP analog channel indices (keys) with their
            corresponding shorthand descriptions (values). (see force_channels)
    """

    # Conditions based on setting the 4 LED bits [first_cue GO_cue], i.e.:
    # (x x x x) = (bottom right, top left, top right, bottom left) LED
    #
    # 1 0 1 0    |  1 0 0 1     = 169
    # right LEDs |  bottom LEDs
    #
    # Note that this is confusing, since the bits are left-right flipped
    # compared to standard notation, where bit 1 is right, bit 4 is left!
    # E.g, right LEDs = bits 1 and 3, this would normally be 0101, however
    # in creating the code we write the bits from left to right: 1010.
    # Then, in calculating the code 169, we again count the rightmost bit
    # as 2^0, the second rightmost as 2^1,...

    # Create a dictionary of trial type codes and corresponding strings
    __grip = {"PG": "1010", "SG": "0101"}
    __force = {"HF": "0110", "LF": "1001"}
    __gf_types = dict(__grip, **__force)

    __grip_first_str = [''.join(i) for i in list(itertools.product(
        list(__grip.keys()), list(__force.keys())))]
#     __grip_only_str = [[''.join(i) for i in list(itertools.product(
#         [list(__grip.keys())[j]], [list(__grip.keys())[j]]))][0] for j in range(2)]
    __grip_only_str = ['PGPG','SGSG']
    __force_first_str = [''.join(i) for i in list(itertools.product(
        list(__force.keys()), list(__grip.keys())))]
#     __force_only_str = [[''.join(i) for i in list(itertools.product(
#         [list(__force.keys())[j]], [list(__force.keys())[j]]))][0] for j in range(2)]
    __force_only_str = ['LFLF','HFHF']

    __list_trial_types = __grip_first_str + __grip_only_str + \
        __force_first_str + __force_only_str

    trial_type_str = {}
    trial_type_codes = {}
    for tt in __list_trial_types:
        trial_type_str[int(__gf_types[tt[:2]] + __gf_types[tt[-2:]], 2)] = tt
        trial_type_codes[tt] = int(__gf_types[tt[:2]] + __gf_types[tt[-2:]], 2)
    trial_type_str[0] = 'n.d.'
    trial_type_codes['n.d.'] = 0
    del tt

    # Create a dictionary of conditions (trial type combinations)
    condition_trial_types_str = {0: [],
                                 1: ['SGHF', 'SGLF', 'PGHF', 'PGLF'],
                                 2: ['HFSG', 'HFPG', 'LFSG', 'LFPG'],
                                 11: ['SGLF', 'PGLF'],
                                 12: ['SGHF', 'PGHF'],
                                 13: ['SGHF', 'SGLF'],
                                 14: ['PGHF', 'PGLF'],
                                 21: ['LFSG', 'LFPG'],
                                 22: ['HFSG', 'HFPG'],
                                 23: ['HFSG', 'LFSG'],
                                 24: ['HFPG', 'LFPG'],
                                 131: ['SGLF'],
                                 132: ['SGHF'],
                                 133: ['SGSG'],
                                 141: ['PGLF'],
                                 142: ['PGHF'],
                                 144: ['PGPG'],
                                 211: ['LFLF'],
                                 213: ['LFSG'],
                                 214: ['LFPG'],
                                 222: ['HFHF'],
                                 223: ['HFSG'],
                                 224: ['HFPG']}

    condition_trial_types_codes = {}
    for c in condition_trial_types_str.keys():
        condition_trial_types_codes[c] = [
            trial_type_codes[i] for i in condition_trial_types_str[c]]
    del i, c

    # Create dictionary to map non-LFP analog channels to channel ID
    force_channels = [137, 138, 139, 140, 141]

    __force_channels_strings = [
        'FSR_PG1', 'FSR_SG1', 'FSR_PG2', 'FSR_SG2', 'LF']

    force_channels_str = {}
    force_channels_idx = {}

    for (i, fc) in enumerate(force_channels):
        force_channels_str[fc] = __force_channels_strings[i]
        force_channels_idx[__force_channels_strings[i]] = fc
    del i, fc

    # Indices of trial id, trial type code and trial events in trial_data
    trial_data_idx = {'trial_id': 0,
                      'trial_id_depr': 1,
                      'performance_code': 2,
                      'trial_type': 3,
                      'det_force': 4,
                      'TS': 5,
                      'FP-ON': 6,
                      'CUE-ON': 7,
                      'CUE-OFF': 8,
                      'GO-ON': 9,
                      'GO-OFF': 10,
                      'SR': 11,
                      'RW': 12,
                      'OT': 13,
                      'OR': 14}

    trial_data_str = {}
    for (k, v) in trial_data_idx.items():
        trial_data_str[v] = k
    del k, v

    # Create a dictionary of event codes, corresponding short-hands and
    # description
    trial_events = [
        'TS', 'FP-ON', 'CUE-ON', 'CUE-OFF', 'GO-ON', 'GO-OFF',
        'SR', 'RW', 'OT', 'OR']

    trial_events_str = {}
    trial_events_codes = {}

    for i in trial_events:
        trial_events_str[trial_data_idx[i]] = i
        trial_events_codes[i] = trial_data_idx[i]
    del i

    # Create diictionaries of performance_code in trial_data
    performance_code_idx = {'correct_trial': 127,
                            'SR<FP-ON': 64,
                            'SR<CUE-ON': 96,
                            'SR<CUE-OFF': 112,
                            'SR<GO': 120,
                            'SR_missing': 124,
                            'grip_error': 126}

    performance_code_str = {}
    for (k, v) in performance_code_idx.items():
        performance_code_str[v] = k
    del k, v

    # This is the length of the filter applied to the force signals for OT
    # detection
    __OT_OR_filterlen = 10
    # This is the threshold on the derivative of the filtered force signals to
    # detect object touch
    __OT_threshold = 25.0
    # This is the threshold on the derivative of the filtered force signals to
    # detect object release
    __OR_threshold = 2.00

    def __init__(
            self, filename, sorting_version=None, metadata_dir=None,
            print_diagnostic=False):
        """
        Constructor
        """

        # Remember choice whether to print diagnostic messages or not
        self._print_diagnostic = print_diagnostic

        # If the specified .nev file does not exist, use default action and
        # print a warning
        if sorting_version is not None:
            if os.path.isfile(r'' + filename + sorting_version + ".nev") \
                    and os.path.isfile(
                        r'' + filename + sorting_version + "-test.txt"):
                txtpostfix = sorting_version + '-test'
            elif os.path.isfile(r'' + filename + sorting_version + ".nev") \
                    and os.path.isfile(
                        r'' + filename + sorting_version + ".txt"):
                txtpostfix = sorting_version
            else:
                self._diagnostic_print(
                    "Requested file %s .nev not found -- "
                    "auto-selecting best available nev." %
                    (filename + sorting_version))
                sorting_version = None

        # If necessary find latest nev file
        if sorting_version is None:
            txtversions = glob.glob(filename + '*.txt')

            # Is there any txt file?
            if len(txtversions) != 0:
                highest_ver = 0

                # If only one txt file exists, use it. If several files exist,
                # use *-test as best version.
                if len(txtversions) == 1:
                    txtpostfix = txtversions[0].\
                        replace(filename, '').replace('.txt', '')
                else:
                    txtpostfix = '-test'
                    if not os.path.isfile(filename + txtpostfix + ".txt"):
                        raise IOError("Unable to auto-detect sorting version \
                            -- multiple .txt files detected.")

                # Get corresponding nev
                nevversions = glob.glob(filename + txtpostfix + '-*.nev')

                for nevversion_i in nevversions:
                    try:
                        # extract version from filename
                        # (last 2 digits before ".nev")
                        # Ex.: -6 -5 -4 -3 -2 -1
                        #       0  5  .  n  e  v
                        vernum = int(nevversion_i[-6:-4])
                    except ValueError:
                        pass
                    else:
                        # If this version is higher, then save the sorting
                        # version postfix
                        if vernum > highest_ver:
                            sorting_version = \
                                nevversion_i.replace(
                                    filename, '').replace(".nev", "")
                            highest_ver = vernum

                # If a nev file without sorting version is found
                if sorting_version is None and \
                        os.path.isfile(filename + txtpostfix + ".nev"):
                    sorting_version = txtpostfix
            else:
                txtpostfix = None

            # If still no sorting was found (either no txt file, or a
            # *-test.txt file exists, but no *-test.nev file exists), then use
            # original nev file
            if sorting_version is None:
                if os.path.isfile(filename + ".nev"):
                    sorting_version = ""
                else:
                    raise IOError("Cannot find a nev file for session " +
                                  "%s." % filename)

        # Output which file is used in the end
        self._diagnostic_print("Using nev file: %s.nev" % (filename +
            sorting_version))

        # Initialize file
        BlackrockIO.__init__(self, filename=filename,
                             nev_override=filename + sorting_version,
                             print_diagnostic=self._print_diagnostic)

        # Determine path to sorting file (txt) if it exists
        if txtpostfix is None:
            self.txt_fileprefix = None
        else:
            self.txt_fileprefix = self.associated_fileset + txtpostfix

        # Save location of metadata directory
        self.metadata_dir = metadata_dir

        # Get odML_doc if available
        self.odML_avail = False
        self.odML_extension = None
        self.odML_filename = None
        self.odML_path = None
        self.odML_doc = None
        if metadata_dir is not None:
            # Make sure the specified directory is a metadata directory
            # TODO remove odMLbuild once server is consistent

            if True not in [os.path.exists(os.path.join(metadata_dir,
                                                        'odMLbuild')),
                            os.path.exists(os.path.join(metadata_dir,
                                                        'odML')),
                            os.path.exists(os.path.join(metadata_dir,
                                                        'build')),
                            os.path.exists(os.path.join(metadata_dir,
                                                        'odMLfiles'))]:
                if not os.path.exists(os.path.join(metadata_dir, 'source')):
                    raise IOError(
                        "Invalid metadata directory specified "
                        "(no build and source directories detected)")

            # Build anticipated location of odML file
            odML_filename = os.path.basename(self.filename)
            if os.path.exists(os.path.join(metadata_dir, 'odMLbuild')):
                odML_path = os.path.join(
                    metadata_dir, 'odMLbuild', odML_filename)
            else:
                odML_path = os.path.join(
                    metadata_dir, 'odMLfiles', odML_filename)
            if os.path.exists(odML_path + '.odml'):
                # Default: odML is in metadata_dir/build
                self.odML_avail = True
                self.odML_extension = 'odml'
                self.odML_filename = odML_filename + '.odml'
                self.odML_path = odML_path + '.odml'
                self.odML_doc = odml.tools.xmlparser.load(self.odML_path)
            elif os.path.exists(odML_path + '.json'):
                self.odML_avail = True
                self.odML_extension = 'json'
                self.odML_filename = odML_filename + '.json'
                self.odML_path = odML_path + '.json'
                self.odML_doc = \
                    odml.tools.jsonparser.json(open(self.odML_path))
            elif os.path.exists(metadata_dir + '.odml'):
                # Alternative: odML is directly in metadata_dir
                # TODO: Remove this -- this code cannot be reached anyhow
                self.odML_avail = True
                self.odML_extension = 'odml'
                self.odML_filename = odML_filename + '.odml'
                self.odML_path = metadata_dir + '.odml'
                self.odML_doc = odml.tools.xmlparser.load(self.odML_path)
            elif os.path.exists(metadata_dir + '.json'):
                self.odML_avail = True
                self.odML_extension = 'json'
                self.odML_filename = odML_filename + '.json'
                self.odML_path = metadata_dir + '.json'
                self.odML_doc = \
                    odml.tools.jsonparser.json(open(self.odML_path))
            else:
                self._diagnostic_print(
                    "No matching odML file found in %s. Attempting to "
                    "assemble minimal metadata from source." % metadata_dir)

        # Interpret file
        self.__get_trial_data()
#        self.__calc_trigger_times()
        self.__load_suamua()

    def __extract_object_touch(
            self, analog_channel_1, analog_channel_2, switch_release_time,
            reward_time):
        """
        Calculates the time stamp of object touch based on the thresholding of
        two signals (which are SG or PG force signals)

        Args:
            analog_channel_1 (neo.AnalogSignal):
                First analog signal.
            analog_channel_2 (neo.AnalogSignal):
                Second analog signal.
            switch_release_time (quantities.Quantity):
                Time of switch release.
            reward_time (quantities.Quantity):
                Time of reward.

        Returns:
            integer:
                Time index at which object touch is detected. If no OT was
                detected, None is returned.
        """

        # TODO: Clean up if all is well
#         a1 = analog_channel_1[int((go_time -
#             2000 * pq.ms).rescale(analog_channel_1.sampling_period).base):
#             int((go_time + 2500 * pq.ms).rescale(
#             analog_channel_1.sampling_period).base)].magnitude
#         a1 = scipy.signal.filtfilt(np.ones((self.__OT_OR_filterlen)) /
#             self.__OT_OR_filterlen, np.array([1]), a1)
#         diff_sig1 = np.diff(a1)
#         a2 = analog_channel_2[int((go_time -
#             2000 * pq.ms).rescale(analog_channel_2.sampling_period).base):
#             int((go_time + 2500 * pq.ms).rescale(
#             analog_channel_2.sampling_period).base)].magnitude
#         a2 = scipy.signal.filtfilt(np.ones((self.__OT_OR_filterlen)) /
#             self.__OT_OR_filterlen, np.array([1]), a2)
#         diff_sig2 = np.diff(a2)

#         diff_sig1 = np.diff(analog_channel_1[int((go_time -
#             2000 * pq.ms).rescale(analog_channel_1.sampling_period).base):
#             int((go_time + 2500 * pq.ms).rescale(
#             analog_channel_1.sampling_period).base)].magnitude)
#
#         diff_sig2 = np.diff(analog_channel_2[int((go_time -
#             2000 * pq.ms).rescale(analog_channel_2.sampling_period).base):
#             int((go_time + 2500 * pq.ms).rescale(
#             analog_channel_2.sampling_period).base)].magnitude)

        diff_sig1 = np.diff(analog_channel_1[int((switch_release_time +
            50 * pq.ms).rescale(analog_channel_1.sampling_period).base):
            int(reward_time.rescale(analog_channel_1.sampling_period).base)]
            .magnitude)

        diff_sig2 = np.diff(analog_channel_2[int((switch_release_time +
            50 * pq.ms).rescale(analog_channel_2.sampling_period).base):
            int(reward_time.rescale(analog_channel_2.sampling_period).base)]
            .magnitude)

        if len(diff_sig1) < 1:
            return None

        th1 = np.max(diff_sig1) / self.__OT_threshold
        th2 = np.max(diff_sig2) / self.__OT_threshold

#         diff_sig1 = np.diff(a1[int((2000 * pq.ms - go_time +
#             switch_release_time +
#             50 * pq.ms).rescale(analog_channel_1.sampling_period).base):
#             int((2000 * pq.ms - go_time + reward_time).rescale(
#             analog_channel_1.sampling_period).base)])
#
#         diff_sig2 = np.diff(a2[int((2000 * pq.ms - go_time +
#             switch_release_time +
#             50 * pq.ms).rescale(analog_channel_2.sampling_period).base):
#             int((2000 * pq.ms - go_time + reward_time).rescale(
#             analog_channel_2.sampling_period).base)])

        ind_sig1 = np.nonzero(diff_sig1 > th1)[0]
        ind_sig2 = np.nonzero(diff_sig2 > th2)[0]

        if len(ind_sig1) > 0:
            t_sig1 = int((switch_release_time + 50 * pq.ms +
                pq.Quantity(ind_sig1[0],
                analog_channel_1.sampling_period).rescale(self.nev_unit)).base)
        else:
            t_sig1 = None
        if len(ind_sig2) > 0:
            t_sig2 = int((switch_release_time + 50 * pq.ms +
                pq.Quantity(ind_sig2[0],
                analog_channel_2.sampling_period).rescale(self.nev_unit)).base)
        else:
            t_sig2 = None

#         plt.plot(diff_sig1, 'g--')
#         plt.plot(diff_sig2, 'm--')
#         plt.hlines(th1, 0, len(diff_sig1), 'g')
#         plt.hlines(th2, 0, len(diff_sig2), 'm')
#         plt.vlines(ind_sig1[0], -200, 200)
#         plt.vlines(ind_sig2[0], -200, 200)

        return max([t_sig1, t_sig2])

    def __extract_object_release(self, analog_channel_1, analog_channel_2,
            reward_time):
        """
        Calculates the time stamp of object release based on the thresholding
        of two signals (which are SG or PG force signals)

        Args:
            analog_channel_1 (neo.AnalogSignal):
                First analog signal.
            analog_channel_2 (neo.AnalogSignal):
                Second analog signal.
            reward_time (quantities.Quantity):
                Time of reward.

        Returns:
            integer:
                Time index at which object release is detected. If no OR was
                detected, None is returned.
        """

        diff_sig1 = np.diff(analog_channel_1[
            int(reward_time.rescale(analog_channel_1.sampling_period).base):
            int((reward_time + 500 * pq.ms).rescale(
                analog_channel_1.sampling_period).base)].magnitude)

        diff_sig2 = np.diff(analog_channel_2[
            int(reward_time.rescale(analog_channel_2.sampling_period).base):
            int((reward_time + 500 * pq.ms).rescale(
                analog_channel_2.sampling_period).base)].magnitude)

        if len(diff_sig1) < 1:
            return None

        th1 = np.min(diff_sig1) / self.__OR_threshold
        th2 = np.min(diff_sig2) / self.__OR_threshold

        ind_sig1 = np.nonzero(diff_sig1 < th1)[0]
        ind_sig2 = np.nonzero(diff_sig2 < th2)[0]

        if len(ind_sig1) > 0:
            t_sig1 = int((reward_time +
                pq.Quantity(ind_sig1[0],
                analog_channel_1.sampling_period).rescale(self.nev_unit)).base)
        else:
            t_sig1 = None
        if len(ind_sig2) > 0:
            t_sig2 = int((reward_time +
                pq.Quantity(ind_sig2[0],
                analog_channel_2.sampling_period).rescale(self.nev_unit)).base)
        else:
            t_sig2 = None

#         import matplotlib.pyplot as plt
#         plt.plot(analog_channel_1[
#             int(reward_time.rescale(analog_channel_1.sampling_period).base):
#             int((reward_time + 500 * pq.ms).rescale(
#             analog_channel_1.sampling_period).base)].magnitude / 100.0, 'r')
#
#         plt.plot(analog_channel_2[
#             int(reward_time.rescale(analog_channel_2.sampling_period).base):
#             int((reward_time + 500 * pq.ms).rescale(
#             analog_channel_2.sampling_period).base)].magnitude / 100.0, 'b')
#         plt.plot(diff_sig1, 'r--')
#         plt.plot(diff_sig2, 'b--')
#         plt.hlines(th1, 0, len(diff_sig1), 'r')
#         plt.hlines(th2, 0, len(diff_sig2), 'b')
#         plt.vlines(ind_sig1[0], -200, 200)
#         plt.vlines(ind_sig2[0], -200, 200)
#         plt.show()

        return max([t_sig1, t_sig2])

    def __extract_force(self, analog_channel, delay_on_time, go_on_time):
        """
        Calculates whether a trial is high force (HF) or low force (LF) based
        on the load force signal.

        Args:
            analog_channel (neo.AnalogSignal):
                Load force signal.
            delay_on_time (quantities.Quantity):
                Time of delay cue on.
            go_on_time (quantities.Quantity):
                Time of GO cue on.

        Returns:
            int:
                Number that depends on the detected force type:
                1: LF
                2: HF
        """

        trial_load_force = analog_channel[
            int(delay_on_time.rescale(
                analog_channel.sampling_period).base):
            int(go_on_time.rescale(
                analog_channel.sampling_period).base)].magnitude

#         plt.plot(trial_load_force, 'r')
#         plt.show()

        return 2 if np.median(trial_load_force) > 200.0 else 1

    def __load_force_signals(self):
        """
        """
        # Load analog signals for object touch detection
        if self.nsx_avail != []:
            # find nsx with highest sampling rate with force information
            mnsx = None
            for nsx_i in self.nsx_avail:
                if set(self.force_channels).issubset(
                        set(self.channel_id_nsx[nsx_i])):
                    mnsx = nsx_i

            if mnsx is not None:
                # get signals and filter (smooth) them
                analog_channels_sg1 = super(ReachGraspIO, self).read_block(
                    channel_list=[self.force_channels_idx['FSR_SG1']],
                    nsx=mnsx, units=None, events=False,
                    waveforms=False).segments[0].analogsignals[0]
                analog_channels_sg1 = \
                    analog_channels_sg1.duplicate_with_new_array(
                        scipy.signal.filtfilt(np.ones((self.__OT_OR_filterlen)) /
                    self.__OT_OR_filterlen, np.array([1]),
                    analog_channels_sg1.magnitude))

                analog_channels_sg2 = super(ReachGraspIO, self).read_block(
                    channel_list=[self.force_channels_idx['FSR_SG2']],
                    nsx=mnsx, units=None, events=False,
                    waveforms=False).segments[0].analogsignals[0]
                analog_channels_sg2 = \
                    analog_channels_sg2.duplicate_with_new_array(
                        scipy.signal.filtfilt(np.ones((self.__OT_OR_filterlen)) /
                    self.__OT_OR_filterlen, np.array([1]),
                    analog_channels_sg2.magnitude))

                analog_channels_pg1 = super(ReachGraspIO, self).read_block(
                    channel_list=[self.force_channels_idx['FSR_PG1']],
                    nsx=mnsx, units=None, events=False,
                    waveforms=False).segments[0].analogsignals[0]
                analog_channels_pg1 = \
                    analog_channels_pg1.duplicate_with_new_array(
                        scipy.signal.filtfilt(np.ones((self.__OT_OR_filterlen)) /
                    self.__OT_OR_filterlen, np.array([1]),
                    analog_channels_pg1.magnitude))

                analog_channels_pg2 = super(ReachGraspIO, self).read_block(
                    channel_list=[self.force_channels_idx['FSR_PG2']],
                    nsx=mnsx, units=None, events=False,
                    waveforms=False).segments[0].analogsignals[0]
                analog_channels_pg2 = \
                    analog_channels_pg2.duplicate_with_new_array(
                        scipy.signal.filtfilt(np.ones((self.__OT_OR_filterlen)) /
                    self.__OT_OR_filterlen, np.array([1]),
                    analog_channels_pg2.magnitude))

                analog_channels_lf = super(ReachGraspIO, self).read_block(
                    channel_list=[self.force_channels_idx['LF']],
                    nsx=mnsx, units=None, events=False,
                    waveforms=False).segments[0].analogsignals[0]
                analog_channels_lf = \
                    analog_channels_lf.duplicate_with_new_array(
                        scipy.signal.filtfilt(np.ones((self.__OT_OR_filterlen)) /
                    self.__OT_OR_filterlen, np.array([1]),
                    analog_channels_lf.magnitude))

                # The DP channel is sometimes 142, sometimes 143. Here we check
                # which one it is for this session
                for dpc in [142, 143]:
                    if dpc in self.channel_id_nsx[nsx_i]:
                        self.force_channels.append(dpc)
                        self.__force_channels_strings.append('DP')
                        self.force_channels_idx[
                            self.__force_channels_strings[-1]] = dpc
                        self.force_channels_str[dpc] = \
                            self.__force_channels_strings[-1]
            else:
                analog_channels_sg1 = None
                analog_channels_sg2 = None
                analog_channels_pg1 = None
                analog_channels_pg2 = None
                analog_channels_lf = None
                # In case needed at a later stage:
#                 analog_channels_dp = super(ReachGraspIO, self).read_block(
#                     channel_list=[self.force_channels_idx['DP']],
#                     nsx=mnsx, units=None, events=False,
#                     waveforms=False).segments[0].analogsignals[0]
        else:
            analog_channels_sg1 = None
            analog_channels_sg2 = None
            analog_channels_pg1 = None
            analog_channels_pg2 = None
            analog_channels_lf = None

        return (analog_channels_sg1, analog_channels_sg2, analog_channels_pg1,
                analog_channels_pg2, analog_channels_lf)

    def __get_OT_OR_force(self, fsig_sg1, fsig_sg2, fsig_pg1, fsig_pg2, fsig,
                          pqCUEon, pqGOon, pqSR, pqRW, trial_type):
        """
        """
        # Manually detect OT/OR events
        if self.is_trialtype_class(trial_type, 'SG'):
            OT = self.__extract_object_touch(fsig_sg1, fsig_sg2, pqSR, pqRW)
            OR = self.__extract_object_release(fsig_sg1, fsig_sg2, pqRW)
        elif self.is_trialtype_class(trial_type, 'PG'):
            OT = self.__extract_object_touch(fsig_pg1, fsig_pg2, pqSR, pqRW)
            OR = self.__extract_object_release(fsig_pg1, fsig_pg2, pqRW)
        else:
            OT = 0
            OR = 0

        if OT > 0:
            det_force = self.__extract_force(fsig, pqCUEon, pqGOon)
        else:
            if self.is_trialtype_class(trial_type, 'HF'):
                det_force = 2
            elif self.is_trialtype_class(trial_type, 'LF'):
                det_force = 1
            else:
                det_force = 0

        return OT, OR, det_force

    def __get_digital_triggers(self):
        """
        """
        # Make sure digital marker index list is sorted
        digi_index = np.sort(self._digital_marker_index)

        # Extract only lower bits of digital markers and corresponding times
        trigger_ids = np.bitwise_and(
            self._event_digital_marker[digi_index], 0x00ff)
        trigger_times = self._event_timestamps[digi_index]

        # Now save the time points of ON and OFF triggers for all 8 digits
        on_indices = []
#        on_times = []
        off_indices = []
#        off_times = []
        # digits:
        #     digit 0: LED bottom right
        #     digit 1: LED top left
        #     digit 2: LED top right
        #     digit 3: LED bottom left
        #     digit 4: Trial start (both home pads pressed)
        #     digit 5: Switch release of working hand
        #     digit 6: Fixation point
        #     digit 7: Reward
        for digit_i in xrange(8):
            # Construct a vector that is 1 if a bit goes on, and -1 if a bit
            # goes off. Note that we purposely do not make assumptions about
            # the state of bits at the beginning and end of the vector, i.e.,
            # we do not know if it a bit is on or off at the beginning -- this
            # information is lost. Therefore, it may be that the vector starts
            # with an off event rather than an on event.
            bit_change_vector = np.diff(
                np.bitwise_and(trigger_ids, 2 ** digit_i) / (2 ** digit_i))

            on_indices.append(np.flatnonzero(bit_change_vector == 1) + 1)
#            on_times.append(trigger_times[on_indices[digit_i]])

            off_indices.append(np.flatnonzero(bit_change_vector == -1) + 1)
#            off_times.append(trigger_times[off_indices[digit_i]])

        return trigger_times, on_indices, off_indices

    def __detect_next_eventidx(self, evX_indices, TS_indices, trial_id,
                               trigger_times, refidx=None):
        """
        """
        # get number of det trials
        no_det_trials = len(TS_indices)

        # define refidx via ref_indices if it is not given
        refidx = refidx if refidx else TS_indices[trial_id - 1]

        # indentify all onsets of evXidx larger than ref_eventidx
        list_indices = np.flatnonzero(evX_indices > refidx)

        # if there are any ...
        if len(list_indices) > 0:
            last_trial = False
            # find next evXidx after ref_eventidx
            next_evXidx = evX_indices[list_indices[0]]

            # if not in the last trial ...
            if trial_id < no_det_trials:
                # check if next_evXidx still belongs tosame trial
                next_TSidx = TS_indices[trial_id]
                if (trigger_times[next_evXidx] >= trigger_times[next_TSidx]):
                    # if no state next_evXidx as undetetable
                    next_evXidx = np.inf

        # if there are none ...
        else:
            last_trial = True
            next_evXidx = np.inf

        return next_evXidx, last_trial

    def __detect_CUEidx(self, on_indices, off_indices, trial_id, trigger_times,
                        refidx_on_ind, refidx_off_ind):
        """
        """
        # Find onsets and offsets of the 4 cue LEDs after the trial_start
        LEDon_idx = [sys.maxint] * 4
        LEDoff_idx = [sys.maxint] * 4
        for LEDi in xrange(4):
            # LED led_i is turned on
            LEDon_idx[LEDi], last_trial = \
                self.__detect_next_eventidx(evX_indices=on_indices[LEDi],
                                            TS_indices=on_indices[4],
                                            trial_id=trial_id,
                                            trigger_times=trigger_times,
                                            refidx=refidx_on_ind)
            # LED led_i is turned off
            LEDoff_idx[LEDi], last_trial = \
                self.__detect_next_eventidx(evX_indices=off_indices[LEDi],
                                            TS_indices=on_indices[4],
                                            trial_id=trial_id,
                                            trigger_times=trigger_times,
                                            refidx=refidx_off_ind)

        # Find next onsets of (delay) cue using on- and offsets of the LEDs
        cueon_idx = min(LEDon_idx)
        cueoff_idx = min(LEDoff_idx) if not np.isinf(cueon_idx) else np.inf

        # Determine cue type (grip or force type for trial type)
        if len(np.unique(LEDon_idx)) == 1:
            cueon_idx = np.inf
            cueoff_idx = np.inf
        if not np.isnan(cueon_idx):
            cue_tmp = sorted(enumerate(LEDon_idx), key=lambda x: x[1])
            cuetype = [0] * 4
            cuetype[cue_tmp[0][0]] = 1
            cuetype[cue_tmp[1][0]] = 1
        else:
            cuetype = [np.inf]

        return cueon_idx, cueoff_idx, cuetype, last_trial

    def __get_trial_data(self):
        """
        Calculates the times of different behavioral events from the digital
        and analog triggers, and the analog signal.

        Args:
            None.

        Returns:
            None.
        """
        # get digital marker indices and times
        trigger_times, on_indices, off_indices = self.__get_digital_triggers()

        # get analog force signals
        force_sig_available = True
        fsig_sg1, fsig_sg2, fsig_pg1, fsig_pg2, fsig = \
            self.__load_force_signals()
        if fsig_sg1 is None:
            force_sig_available = False

        # Extract all trials, successful and errors into trialdata.
        # Trial ids:
        #     trial_id:         ids for each detected trial_start,
        #                       error trials are included
        #     trial_ids_dep:    ids for each detected reward,
        #                       error trials are omitted
        # Digital trial events:
        #     TS:         on_indices idx 4
        #     FP-ON:      on_indices idx 6
        #     CUE-ON:     on_indices idx 0-3
        #     CUE-OFF:    off_indices idx 0-3
        #     GO-ON:      on_indices idx 0-3
        #     GO-OFF:     off_indices idx 0-3
        #     SR:         on_indices idx 5
        #     RW:         on_indices idx 7
        # Analog trial events (in nsx file):
        #     OT:    channel_id 137-140
        #     OR:    channel_id 137-140
        # Performance code:
        #     127:    successful trial
        #     64:     HPR before FP onset
        #     96:     HPR before cue onset
        #     112:    HPR before cue offset
        #     120:    HPR before GO
        #     124:    No HPR after GO
        #     126:    No Reward, Grip Error

        # Detect number of TRIAL STARTS:
        no_detTS_trials = len(on_indices[4])

        # Create trialdata array:
        trialdata = np.zeros((no_detTS_trials,
                              len(self.trial_data_idx.keys())), dtype=int)

        if len(trialdata) > 0:
            # start counting for old trial ids (only correct trials)
            id_correct_trial = 0
            # test if first trial is without trial start
            first_TS = on_indices[4][0]
            if len(on_indices[7]) > 0:
                first_RW = on_indices[7][0]
                if (trigger_times[first_TS] >= trigger_times[first_RW]):
                    id_correct_trial += 1

        # use each trial start to detect remaining trialdata
        for (i, TS_idx) in enumerate(on_indices[4]):
            # Find next onsets of fixation (6) after the trial_start
            FPon_idx, last_trial = \
                self.__detect_next_eventidx(evX_indices=on_indices[6],
                                            TS_indices=on_indices[4],
                                            trial_id=i + 1,
                                            trigger_times=trigger_times)
            if last_trial:
                self._diagnostic_print("Trial out of rec time - omitting!")
                trialdata = trialdata[:-1]
                continue

            # Find onsets and offsets of the 4 cue LEDs after the trial_start
            if not np.isinf(FPon_idx):
                CUEon_idx, CUEoff_idx, delay_cue, last_trial = \
                    self.__detect_CUEidx(on_indices=on_indices,
                                         off_indices=off_indices,
                                         trial_id=i + 1,
                                         trigger_times=trigger_times,
                                         refidx_on_ind=FPon_idx,
                                         refidx_off_ind=FPon_idx)
                if last_trial:
                    self._diagnostic_print("Trial out of rec time - omitting!")
                    trialdata = trialdata[:-1]
                    continue
            else:
                CUEon_idx = np.inf
                CUEoff_idx = np.inf
                delay_cue = [np.inf]

            # if a first cue was found, find the second cue
            if not np.isinf(CUEon_idx) and not np.isinf(CUEoff_idx):
                GOon_idx, GOoff_idx, go_cue, last_trial = \
                    self.__detect_CUEidx(on_indices=on_indices,
                                         off_indices=off_indices,
                                         trial_id=i + 1,
                                         trigger_times=trigger_times,
                                         refidx_on_ind=CUEoff_idx,
                                         refidx_off_ind=CUEoff_idx)
                if last_trial:
                    self._diagnostic_print("Trial out of rec time - omitting!")
                    trialdata = trialdata[:-1]
                    continue

                if not np.isinf(GOon_idx) and \
                        trigger_times[GOon_idx] < trigger_times[CUEoff_idx]:
                    GOon_idx = np.inf
                    GOoff_idx = np.inf
                    CUEoff_idx = np.inf

            # if no first cue was found define second cue also as undetectable
            else:
                GOon_idx = np.inf
                GOoff_idx = np.inf
                go_cue = [np.inf]

            # Find switch release (5) after trial start
            if not np.isinf(GOon_idx):
                SR_idx, last_trial = \
                    self.__detect_next_eventidx(evX_indices=on_indices[5],
                                                TS_indices=on_indices[4],
                                                trial_id=i + 1,
                                                trigger_times=trigger_times)
                if last_trial:
                    self._diagnostic_print("Trial out of rec time - omitting!")
                    trialdata = trialdata[:-1]
                    continue
            else:
                SR_idx = np.inf

            # Find next reward (7) after the trial start
            RW_idx, last_trial = \
                self.__detect_next_eventidx(evX_indices=on_indices[7],
                                            TS_indices=on_indices[4],
                                            trial_id=i + 1,
                                            trigger_times=trigger_times)
            if last_trial:
                self._diagnostic_print("Trial out of rec time - omitting!")
                trialdata = trialdata[:-1]
                continue

            # Calculate trial type. Note that the bits of the go and delay
            # cue are each flipped left to right! E.g., bit 1 of the go cue
            # counts as 2^3, bit 2 as 2^2, bit 1 as 2^1, and bit 0 as 2^0
            trial_type = 0
            if not np.isinf(GOon_idx) and not np.isinf(CUEon_idx):
                for a, b in enumerate(go_cue):
                    trial_type += b * 2 ** (3 - a)
                for a, b in enumerate(delay_cue):
                    trial_type += b * 2 ** ((3 - a) + 4)

            # Add every event and trial information to trialdata_all list
            TS = trigger_times[TS_idx]
            FPon = 0 if np.isinf(FPon_idx) else trigger_times[FPon_idx]
            CUEon = 0 if np.isinf(CUEon_idx) else trigger_times[CUEon_idx]
            CUEoff = 0 if np.isinf(CUEoff_idx) else trigger_times[CUEoff_idx]
            GOon = 0 if np.isinf(GOon_idx) else trigger_times[GOon_idx]
            GOoff = 0 if np.isinf(GOoff_idx) else trigger_times[GOoff_idx]
            SR = 0 if np.isinf(SR_idx) else trigger_times[SR_idx]
            RW = 0 if np.isinf(RW_idx) else trigger_times[RW_idx]

            # get force and time of object touch and object release
            if 0 not in [TS, CUEon, GOon, SR, RW] and force_sig_available:
                OT, OR, det_force = self.__get_OT_OR_force(
                    fsig_sg1, fsig_sg2, fsig_pg1, fsig_pg2, fsig,
                    pq.Quantity(CUEon, self.nev_unit),
                    pq.Quantity(GOon, self.nev_unit),
                    pq.Quantity(SR, self.nev_unit),
                    pq.Quantity(RW, self.nev_unit),
                    trial_type)
            else:
                OT, OR, det_force = 0, 0, 0

            if OT is None:
                self._diagnostic_print("Trial with SR and RW in quick" +
                                       " succession - omitting!")
                trialdata = trialdata[:-1]
                id_correct_trial += 1
                continue

            # Error classification
            ev_binary = ''
            for ev in [TS, FPon, CUEon, CUEoff, GOon, SR, RW]:
                ev_binary += '1' if ev > 0 else '0'
            performance_code = int(ev_binary, 2)

            # define trial_id
            trial_id = i + 1

            # define trial_id_depr
            if performance_code == 127:
                id_correct_trial += 1
                trial_id_depr = id_correct_trial
            else:
                trial_id_depr = -1

            for tdi, td in enumerate([trial_id, trial_id_depr,
                                      performance_code, trial_type,
                                      det_force, TS, FPon, CUEon, CUEoff, GOon,
                                      GOoff, SR, RW, OT, OR]):
                trialdata[trial_id - 1][tdi] = td

            tdprint = []
            td = trialdata[trial_id - 1]
            for it in sorted(self.trial_data_idx.items(), key=lambda x: x[1]):
                tdname = "%s: " % it[0]
                if it[0] == 'trial_type':
                    tddata = "%s    " % self.trial_type_str[td[it[1]]]
                elif it[0] == 'det_force':
                    tddata = "%s    " % ['n.d.', 'LF', 'HF'][td[it[1]]]
                else:
                    tddata = "%s    " % str(td[it[1]])
                tdprint.append(tdname + tddata)

            self._diagnostic_print(''.join(tdprint))

            if det_force == 2 and self.is_trialtype_class(trial_type, 'LF'):
                self._diagnostic_print(
                    '   -> LF type trial with high load force detected!')
            if det_force == 1 and self.is_trialtype_class(trial_type, 'HF'):
                self._diagnostic_print(
                    '   -> HF type trial with low load force detected!')

        self.trial_data = trialdata

        # determine condition and trial numbers if trial_data were generated
        self.condition = 999
        if len(self.trial_data) > 0:
            # trial type index in trial_data
            tt_idx = self.trial_data_idx['trial_type']

            # trial types occurring in one cue task
            one_cue_tt = [85, 102, 153, 170]

            # Is this a session with one cue or two cues?
            self.is_one_cue_task = bool(np.unique(
                [i in one_cue_tt for i in self.trial_data[:, tt_idx]])[0])
            self.is_two_cues_task = True if not self.is_one_cue_task \
                else False

            # Determine all occurring trial types
            self.trial_types = np.sort(
                np.unique(self.trial_data[:, tt_idx].T)).tolist()

            # Determine condition
            self.condition = 0
            for cond_i in self.condition_trial_types_codes.keys():
                if self.trial_types == self.condition_trial_types_codes[
                        cond_i]:
                    self.condition = cond_i

            # Save total number of trials
            self.num_total_trials = len(self.trial_data)
            # Save number of correct and grip error trials
            pc_idx = self.trial_data_idx['performance_code']
            # number of correct trials
            pc_ct = self.performance_code_idx['correct_trial']
            self.num_correct_trials = \
                len(np.nonzero(self.trial_data[:, pc_idx] == pc_ct)[0])
            # number of grip error trials
            pc_ge = self.performance_code_idx['grip_error']
            self.num_grip_error_trials = \
                len(np.nonzero(self.trial_data[:, pc_idx] == pc_ge)[0])

    def __load_suamua(self):
        """
        This function loads the text file associated with a sorted session and
        extracts the number of SUA, their IDs, and the channel of the MUA.

        Args:
            None.

        Returns:
            None.
        """

        self.__sua_ids = []
        self.__mua_ids = []

        # Text file that contains the mapping from 10X10 grid to electrode ID
        # along the linear index (from lower left to upper right, each row
        # starting from left, channel_i.e., Alexa's mapping):
        # 1: Blackrock ID of channel
        # 2: number of SUA
        # 3: if 0, no MUA; if >0: ID of MUA in nev file
        #    (non-sorted spike waveforms)

        if self.txt_fileprefix != None:
            filename_txt = self.txt_fileprefix + '.txt'

            # Load file
            nid_file_withNAN = np.loadtxt(filename_txt)

            # Remove all NaN entries and sort
            mask = np.logical_not(np.isnan(nid_file_withNAN).any(axis=1))
            nid_file_noNAN = nid_file_withNAN[mask].astype(int)
            nid_file_noNAN = nid_file_noNAN[nid_file_noNAN[:, 0].argsort()]

            for channel_i in xrange(96):
                pos = np.nonzero(nid_file_noNAN[:, 0] == channel_i + 1)

                # Make sure only one line per channel exists in the txt file
                if len(pos) != 1:
                    self.__sua_ids = [[] for _ in xrange(96)]
                    self.__mua_ids = [[] for _ in xrange(96)]

                    self.is_sorted = False

                    self._diagnostic_print("MUA/SUA file %s is corrupt -- "
                        "multiple entries for channel %i present." %
                        (filename_txt, channel_i + 1))

                # Make sure the MUA ID does not overlap with the SUA IDs
                elif nid_file_noNAN[pos, 2] in \
                        range(1, nid_file_noNAN[pos, 1] + 1):
                    self.__sua_ids = [[] for _ in xrange(96)]
                    self.__mua_ids = [[] for _ in xrange(96)]

                    self.is_sorted = False

                    self._diagnostic_print("MUA/SUA file %s is corrupt -- "
                        "MUA and SUA information for channel %i is "
                        "corrupt." % (filename_txt, channel_i + 1))

                # Load MUA and SUA ids
                else:
                    pos = pos[0]
                    self.__sua_ids.append(range(1, nid_file_noNAN[pos, 1] + 1))
                    if nid_file_noNAN[pos, 2] != 0:
                        self.__mua_ids.append(list(nid_file_noNAN[pos, 2]))
                    else:
                        self.__mua_ids.append([])

                    self.is_sorted = True
        else:
            self.__sua_ids = [[] for _ in xrange(96)]
            self.__mua_ids = [[] for _ in xrange(96)]

            self.is_sorted = False

            self._diagnostic_print("No distinction between MUA and SUA - no " +
                                   "txt file found.")

    def get_mua_ids(self, electrode):
        """
        Returns a list of MUA IDs recorded on a specific electrode.

        Args:
            channel_id (int):
                Electrode number 1-96 for which to extract neuron IDs.

        Returns:
            list
                List containing all MUA IDs of the specified electrode.
                If no sorting (i.e., no .txt file) exists, the return is [].
        """
        if electrode < 1 or electrode > 96:
            raise Exception("Invalid electrode ID specified.")

        return self.__mua_ids[electrode - 1]

    def get_sua_ids(self, electrode):
        """
        Returns a list of SUA IDs recorded on a specific electrode.

        Args:
            channel_id (int):
                Electrode number 1-96 for which to extract neuron IDs.

        Returns:
            list
                List containing all SUA IDs of the specified electrode.
                If no sorting (i.e., no .txt file) exists, the return is [].
        """

        if electrode < 1 or electrode > 96:
            raise Exception("Invalid electrode ID specified.")

        return self.__sua_ids[electrode - 1]

    def is_trialtype_class(self, trial_type, req):
        """
        Returns a True if a given trial type code is part of the class of
        SG, PG, LF, or HF trials, respectively.

        Args:
            trial_type (int):
                Supplied trial type code
            req (str):
                String of required trial type to test for. Can be either of
                'SG', 'PG', 'LF', or 'HF'.

        Returns:
            bool:
                True, if trial_type belongs to the class of trials specified by
                req.
        """

        if trial_type in self.trial_type_str.keys() and \
           req in self.trial_type_str[trial_type]:
            return True
        return False

    def get_trial_ids(
            self, trial_types=None, performance_codes=None,
            trial_id_type='trial_id'):
        """
        Returns the trial IDs of those trials matching a specific list of trial
        types.

        Args:
            trial_types (list):
                List of integers specifying the trial code (see class
                documentation) or strings of the form 'LFPG', ... of those
                trials whose ID will be returned. If None, all trial IDs
                defined by trial_id_type are returned.
                Default: None.
            performance_codes (list):
                List of integers specifying the trial performance as code (see
                class documentation) or strings of those trials whose ID will
                be returned, i.e. specify [127] or ['correct_trial'] for
                trials with correct trial performance. If None, all trial IDs
                defined by trial_id_type are returned.
                Default: None.
            trial_id_types (string):
                Defines which type of id numbering is returned:
                    trial_id: determined by detection of trial starts
                              (all error trials included)
                    trial_id_depr: determined by detection of reward
                                   (basically no error trials included)
                Default: trial_id.

        Returns:
            List of trial IDs matching the specified trial types, and trial
            performances.

        Example:
            >>> trials=rgobject.get_trial_ids()
            Returns all trial IDs defined by detected trial starts.

            >>> trials=rgobject.get_trial_ids(trial_types='SGLF',
            >>>                               performance_codes=127)
            Returns trial IDs of correct trials of type 'SGLF'.
        """
        if len(self.trial_data) > 0:

            tid_idx = self.trial_data_idx[trial_id_type]
            idlistall = self.trial_data[:, tid_idx].T
            idlist = list(idlistall)

            if trial_types:
                tt_idx = self.trial_data_idx['trial_type']
                # If only a single string or int is given, convert to a list
                if type(trial_types) != list:
                    trial_types = list([trial_types])

                # find ids with given trial types
                idlisttt = []
                for tt in trial_types:
                    # change tt from string to integer code if necessary
                    if isinstance(tt, str):
                        tt = self.trial_type_codes[tt]
                    # create mask for tt and get corresponding trial ids
                    tt_mask = (self.trial_data[:, tt_idx] == tt)
                    idlisttt.extend(idlistall[tt_mask])
                idlist = idlisttt

            if performance_codes:
                pc_idx = self.trial_data_idx['performance_code']
                # If only a single string or int is given, convert to a list
                if type(performance_codes) != list:
                    performance_codes = list([performance_codes])

                # find ids with given performance codes
                idlistpc = []
                for pc in performance_codes:
                    # change tt from string to integer code if necessary
                    if type(pc) is str:
                        pc = self.performance_code_idx[pc]
                    # create mask for pc and get corresponding trial ids
                    pc_mask = (self.trial_data[:, pc_idx] == pc)
                    idlistpc.extend(idlistall[pc_mask])
                idlist = idlistpc

            if trial_types and performance_codes:
                idlist = sorted(list(set(idlisttt) & set(idlistpc)))

        else:
            idlist = []

        return idlist

    def read_block(self, lazy=False, cascade=True, n_starts=[None],
                   n_stops=[None], channel_list=None, nsx=None, units=None,
                   events=False, waveforms=False, corrections=False):
        """
        Reads file contents as a neo Block.

        The Block contains one Segment for each entry in zip(n_starts,
        n_stops). If these parameters are not specified, the default is
        to store all data in one Segment.

        The Block contains one RecordingChannelGroup per channel.

        Args:
            lazy (bool):
            cascade (bool):
            n_starts (list):
                List of starting times as Quantity objects of each Segment,
                where the beginning of the file is time 0. If a list entry is
                specified as None, the beginning of the file is used. Default:
                [None].
            n_stops (list):
                List of corresponding stop times as Quantity objects of each
                Segment. If a list entry is specified as None, the end of the
                file is used. Default: [None].
            channel_list (list):
                List of channels to consider in loading (Blackrock channel
                IDs). The neural data channels are 1 - 128. The analog inputs
                are 129 - 255. If None is specified, all channels are loaded.
                Default: None.
            nsx (list):
                List of integers X (between 0 and 9) specifying the .nsX file
                to load analog data from channels in channel_list from. If an
                empty list is given, all available nsX files are loaded.If None
                is specified, no analog data is read. Blackrock defined
                sampling frequencies:
                    .ns6 - 30000 Hz (no digital filter)
                    .ns5 - 30000 Hz
                    .ns4 - 10000 Hz
                    .ns3 -  2000 Hz
                    .ns2 -  1000 Hz
                    .ns1 -   500 Hz
                Default: None.
            units (list or dictionary):
                Specifies the unit IDs to load from the data. If an empty list
                is specified, all units of all channels in channel_list are
                loaded. If a list of units is given, all units matching one of
                the IDs in the list are loaded from all channels in
                channel_list. If a dictionary is given, each entry with a key N
                contains a list of unit IDs to load from channel N. If None is
                specified, no units are loaded. Blackrock definition of units:
                        0     unsorted
                        1-16  sorted units
                        255   "noise"
                Default: None.
            events (boolean):
                If True, all digital and analog events are inserted as
                EventArrays.
                Default: False
            waveforms (boolean):
                If True, waveforms of each spike are loaded in the SpikeTrain
                objects. Default: False
            corrections (boolean):
                If True, gap file (corrections.txt) is loaded and known gaps are
                corrected

        Returns:
            neo.Block
                neo Block object containing the data.
                Attributes:
                    name: short session name
                    file_origin: session name
                    rec_datetime: date and time
                    description: string "Blackrock format file"
                Annotations:
                    task_type (str): name of performed task type
                    condition (int): code for performed condition of trial type combination
                    random_grip (bool): if True, trial type order according to grip is random
                    random_force (bool): if True, trial type order according to force is random
                    block_size (int): number of repetitions of same grip and/or same force type, if grip and/or force are not random
                    elid_list (list): electrode ids of blackrock
                    elid_list_ca (list): electrode ids of blackrock ordered
                                         according to connector alignment
                                         (top row (left to right column) to
                                         bottom row)
                    elid_list_ba (list): electrode ids of blackrock ordered
                                         according to brain alignment
                                         (top row (left to right column) to
                                         bottom row)

                The neo Block contains the following neo structures:
                Segment
                    For each pair of n_start and n_stop values, one Segment is
                    inserted.
                    Attributes:
                        name: string of the form "Segment i"
                        file_origin: session name
                        rec_datetime: date and time
                        index: consecutive number of the Segment

                RecordingChannelGroup
                    For each recording channel, one RecordingChannelGroup object
                    is created.
                    Attributes:
                        name: string of the form "Channel i"
                        file_origin: session name
                        rec_datetime: date and time
                        channel_indexes: list of one element, which is the integer channel number

                RecordingChannel
                    For each recording channel, one RecordingChannel object is
                    created as child of the respective RecordingChannelGroup.
                    Attributes:
                        name: string of the form "Channel i"
                        file_origin: session name
                        rec_datetime: date and time
                        channel_indexes: list of one element, which is the integer channel number

                AnalogSignal
                    For each loaded analog signal, one AnalogSignal object is created per Segment.
                    Attributes:
                        name: string of the form "Analog Signal Segment i, Channel j, NSX x"
                        file_origin: session name
                    Annotations:
                        channel_id (int): channel number
                        electrode_id (int): electrode id in blackrock
                        ca_id (int): electrode id when array is connector aligned (1-100)
                        ba_id (int): electrode id when array is brain aligned (1-100)
                        rejLFC (bool): True if electrode needs to be rejected
                                       for low frequency components
                        rejIFC (bool): True if electrode needs to be rejected
                                       for intermediate frequency components
                        rejHFC (bool): True if electrode needs to be rejected
                                       for high frequency components
                        file_nsx (int): number X of .nsX file

                Unit
                    For each unit, one Unit structure is created.
                    Attributes:
                        name: string of the form "Channel j, Unit u"
                        file_origin: session name
                        channel_indexes: list of one element, which is the integer channel number
                    Annotations:
                        channel_id (int): channel number
                        electrode_id (int): electrode id in blackrock
                        ca_id (int): electrode id when array is connector aligned (1-100)
                        ba_id (int): electrode id when array is brain aligned (1-100)
                        rejLFC (bool): True if electrode needs to be rejected
                                       for low frequency components
                        rejIFC (bool): True if electrode needs to be rejected
                                       for intermediate frequency components
                        rejHFC (bool): True if electrode needs to be rejected
                                       for high frequency components
                        unit_id (int): unit number
                        sua (bool): True if unit is a SUA
                        mua (bool): True if unit is a MUA

                SpikeTrain
                    For each Unit and each Segment, one SpikeTrain is created.
                    Waveforms of spikes are inserted into the waveforms property of the
                    respective SpikeTrain objects. Individual Spike objects are not
                    used.
                    Attributes:
                        name: string of the form "Segment i, Channel j, Unit u"
                        file_origin: session name
                        dtype: int (original time stamps are save in units of nev_unit)
                        sampling_rate: Waveform time resolution
                    Annotations:
                        channel_id (int): channel number
                        electrode_id (int): electrode id in blackrock
                        ca_id (int): electrode id when array is connector aligned (1-100)
                        ba_id (int): electrode id when array is brain aligned (1-100)
                        rejLFC (bool): True if electrode needs to be rejected
                                       for low frequency components
                        rejIFC (bool): True if electrode needs to be rejected
                                       for intermediate frequency components
                        rejHFC (bool): True if electrode needs to be rejected
                                       for high frequency components
                        unit_id (int): unit number
                        sua (bool): True if unit is a SUA
                        mua (bool): True if unit is a MUA

                EventArray
                    Per Segment, for digital or analog markers with a common ID,
                    one EventArray is created containing the time stamps of all
                    occurrences of that specific marker.
                    In addition, Reach-grasp specific events are inserted for each
                    event in ReachGraspIO.events_str.
                    Attributes:
                        name: string of the form "Digital Marker m" / "Analog Channel i Marker m"
                        file_origin: session name
                        dtype: int (original time stamps are save in units of nev_unit)
                        labels: contains a string of the ID per time stamp
                    Annotations:
                        marker_ID (int): marker ID
                        digital_marker (bool): True if it is a digital marker
                        analog_marker (bool):  True if it is an analog marker
                        analog_channel (int): For analog markers, determines the respective channel (0-4)

            Notes:
                For Segment and SpikeTrain objects, if t_start is not
                specified, it is assumed 0. If t_stop is not specified, the
                maximum time of all analog samples available, all available
                spike time stamps and all trigger time stamps is assumed.
        """

        # reading correction parameters from 'corrections.txt' file and saving them
        # gap_corrections = [gap_start_bin,gap_size_bins]
        gap_corrections = [None, None]
        if corrections:
            try:
                correction_file = open(
                    os.path.dirname(__file__) + '/corrections.txt', 'r')
                for line in correction_file:
                    if os.path.basename(self.filename) in line:
                        numbers = [int(s) for s in line.split() if s.isdigit()]
                        if len(numbers) == 2:
                            gap_corrections = numbers * \
                                np.array(
                                    1.0) * pq.CompoundUnit('1.0/%i*s' % (self.timestamp_res))
                        else:
                            warnings.warn(
                                'Wrong number of integers in corrections.txt for session %s' % os.path.basename(self.filename))
                        break
                correction_file.close()
            except IOError:
                warnings.warn('No file "corrections.txt" found.')

            # correcting n_starts and n_stops for gap
            # listify if necessary
            n_starts_c = n_starts if type(n_starts) == list else [n_starts]
            n_stops_c = n_stops if type(n_stops) == list else [n_stops]

            # shift start and stop times to allow gap correction if gap is
            # known
            if gap_corrections[0] != None:
                for time_list in [n_starts_c, n_stops_c]:
                    # iterate over all n_start and n_stops
                    for i in range(len(time_list)):
                        if time_list[i] >= gap_corrections[0]:
                            time_list[i] += gap_corrections[1]

        # Load neo block
        block = BlackrockIO.read_block(self, lazy=lazy, cascade=cascade,
                                       n_starts=n_starts, n_stops=n_stops,
                                       channel_list=channel_list, nsx=nsx,
                                       units=units, events=events,
                                       waveforms=waveforms)

        # post correct gaps if gap is known
        if corrections and gap_corrections[0] != None:
            # correct alignment
            for i in range(len(block.segments)):

                #                # use t_start if specified
                #                t_start = None
                #                if n_starts:
                #                    t_start = n_starts[i]
                #
                #                # use t_stop if specified
                #                t_stop = None
                #                if n_stops:
                #                    t_stop = n_stops[i]

                # adjust spiketrains
                for j in range(len(block.segments[i].spiketrains)):
                    st = block.segments[i].spiketrains[j]

                    # adjusting t_start
                    if st.t_start >= gap_corrections[0] + gap_corrections[1]:
                        st.t_start -= gap_corrections[1]

                    # correct for gap
                    st = st - ((st > gap_corrections[0]) * gap_corrections[1])

#                    # discard spikes before t_start
#                    if t_start:
#                        idx_valid = np.nonzero(st >= t_start)[0]
#                        if len(idx_valid):
#                            st = st[idx_valid[0]:]
#
#                    # discard spikes after t_stop
#                    if t_stop:
#                        idx_invalid = np.nonzero(st >= t_stop)[0]
#                        if len(idx_invalid):
#                            st = st[:idx_invalid[0]]

                    # shallow copy from original spiketrain (annotations,
                    # waveforms, etc.)
                    st.__dict__ = block.segments[
                        i].spiketrains[j].__dict__.copy()

                    # adjusting t_stop
                    if st.t_stop >= gap_corrections[0] + gap_corrections[1]:
                        st.t_stop -= gap_corrections[1]

                    # link block to new spiketrain
                    block.segments[i].spiketrains[j] = st

        # TODO: odML <-> nev consistency checks? Low priority
        # condition, trials, SUA IDs, MUA IDs

        # Block annotations
        if self.condition is not None:
            block.annotate(condition=self.condition)

        # Add annotations of odML meta data info
        if self.odML_avail:
            # Annotate electrode id lists
            ff = lambda x: x.name.startswith('Electrode_')
            sobjs = [s for s in self.odML_doc.itersections(filter_func=ff)]

            # Annotate electrode id lists
            elid_list = []  # general available ids
            elid_list_ca = np.zeros((100)) - 1  # for connector alignement
            elid_list_ba = np.zeros((100)) - 1  # for brain alignement
            for s in sobjs:
                elid_list.append(s.properties['ID'].value.data)
                elid_list_ca[s.properties['IDca'].value.data - 1] = \
                    s.properties['ID'].value.data
                elid_list_ba[s.properties['IDba'].value.data - 1] = \
                    s.properties['ID'].value.data
            block.annotate(elid_list=sorted(elid_list))
            block.annotate(elid_list_ca=elid_list_ca)

            # Brain-aligned will be dropped
#             block.annotate(elid_list_ba=elid_list_ba)

            # Annotate performed task type
            ff = lambda x: x.name in 'TaskType' and \
                'Subsession' == x.parent.name
            pobjs = [p for p in self.odML_doc.iterproperties(filter_func=ff)]
            vdata = [p.value.data for p in pobjs if p.value.data]
            block.annotate(task_type=vdata[0])

            # Annotate trial type order (random, block, size of block)
            cue_task_types = ['Observation', 'OneCue', 'TwoCues']
            if block.annotations['task_type'] in cue_task_types:
                ff = lambda x: x.name in ['OrderGrip', 'OrderForce'] and \
                    'TrialTypeSettings' in x.parent.name
                pobjs = [
                    p for p in self.odML_doc.iterproperties(filter_func=ff)]
                block.annotate(
                    random_grip=['random' == p.value.data for p in
                                 pobjs if p.name == 'OrderGrip'][0])
                block.annotate(
                    random_force=['random' == p.value.data for p in
                                  pobjs if p.name == 'OrderForce'][0])

                if False in [block.annotations['random_grip'],
                             block.annotations['random_force']]:
                    ff = lambda x: x.name == 'BlockSize' and \
                        'TrialTypeSettings' in x.parent.name
                    pobjs = [p for p in
                             self.odML_doc.iterproperties(filter_func=ff)]
                    block.annotate(block_size=[p.value.data for p in pobjs][0])

        elif self.metadata_dir is not None:
            # If no odML is available, but metadata directory was specified,
            # read a bare minimum of information

            # Annotate available electrode id list for connector alignement
            einfo = rgodml.metadata_io.load_blackrock_electrodes_info(
                os.path.join(
                    self.metadata_dir, 'source', 'blackrock',
                    'electrodes_info.txt'))
            elid_list = sorted([einfo[k]['ID']['data'] for k in einfo.keys()
                                if einfo[k]['ID']['data'] > 0])
            elid_list_ca = list(np.ones(100, dtype=int) * -1)
            for elid in elid_list:
                elid_str = 'Electrode_%03d' % elid
                elid_list_ca[einfo[elid_str]['IDca']['data'] - 1] = elid
            block.annotate(elid_list=elid_list)
            block.annotate(elid_list_ca=elid_list_ca)
            # Annotate available electrode id list for brain alignement
            # Removed -- brainaligned will be dropped
#             einfo = rgodml.metadata_io.load_brain_aligned_elids(
#                 os.path.join(
#                     self.metadata_dir, 'source', 'monkey',
#                     'brain_aligned_elids.txt'))
#             elid_list_ba = list(np.ones(100, dtype=int) * -1)
#             for elid in elid_list:
#                 elid_str = 'Electrode_%03d' % elid
#                 elid_list_ba[einfo[elid_str]['IDba']['data'] - 1] = elid
#             block.annotate(elid_list_ba=elid_list_ba)

        # Add interpreted events to block
        if events:
            for seg_i, seg in enumerate(block.segments):
                # Find index of reach-grasp event in this segment
                for event_i in self.trial_events_str.keys():
                    if n_starts[seg_i] == None:
                        tstart = 0
                    else:
                        tstart = n_starts[seg_i] / self.nev_unit
                    if n_stops[seg_i] == None:
                        tstop = sys.maxint
                    else:
                        tstop = n_stops[seg_i] / self.nev_unit
                    idx = np.nonzero(np.logical_and(
                        self.trial_data[:, event_i] >= tstart,
                        self.trial_data[:, event_i] < tstop))[0]
                    ev = neo.EventArray(
                        times=pq.Quantity(self.trial_data[idx, event_i],
                            units=self.nev_unit,
                            dtype="int"),
                        labels=np.tile(self.trial_events_str[event_i],
                            (len(idx))),
                        name=self.trial_events_str[event_i],
                        file_origin=self.associated_fileset,
                        marker_id=event_i,
                        digital_marker=True,
                        analog_marker=False,
                        analog_channel=0,
                        event_name=self.trial_events_str[event_i])
                    seg.eventarrays.append(ev)

        # Add annotations of spike sorting and odML meta data info
        for seg in block.segments:
            # Add annotations to analogsignals
            for asig in seg.analogsignals:
                if asig.annotations['channel_id'] <= 100:
                    el_id = asig.annotations['channel_id']
                    asig.annotations['electrode_id'] = el_id
                    if self.odML_avail:
                        # Annotate connector and brain aligned id
                        idx = block.annotations['elid_list'].index(el_id)
                        asig.annotations['ca_id'] = \
                            block.annotations['elid_list_ca'][idx]
                        asig.annotations['ba_id'] = \
                            block.annotations['elid_list_ba'][idx]

                        # Annotate if electrode should be rejected
                        ff = lambda x: x.name == 'RejElectrodes' and \
                            'RejectionsLFP' in x.get_path()
                        pobjs = [p for p in
                                 self.odML_doc.iterproperties(filter_func=ff)]
                        vdata = [v.data for v in pobjs[0].values]
                        fc = pobjs[0].parent.name
                        if el_id in vdata:
                            asig.annotations['rej' + fc] = True
                        else:
                            asig.annotations['rej' + fc] = False
                    elif self.metadata_dir is not None:
                        # Annotate connector and brain aligned id
                        asig.annotations['ca_id'] = block.annotations[
                            'elid_list_ca'].index(el_id) + 1

            # Add annotations to spiketrains
            for st in seg.spiketrains:
                if st.annotations['channel_id'] <= 100:
                    el_id = st.annotations['channel_id']
                    st.annotations['electrode_id'] = el_id
                    if self.odML_avail:
                        # Annotate connector and brain aligned id
                        idx = block.annotations['elid_list'].index(el_id)
                        st.annotations['ca_id'] = \
                            block.annotations['elid_list_ca'][idx]
                        st.annotations['ba_id'] = \
                            block.annotations['elid_list_ba'][idx]

                        # Annotate if electrode should be rejected
                        ff = lambda x: x.name == 'RejElectrodes' and \
                            'RejectionsLFP' in x.get_path()
                        pobjs = [p for p in
                                 self.odML_doc.iterproperties(filter_func=ff)]
                        vdata = [v.data for v in pobjs[0].values]
                        fc = pobjs[0].parent.name
                        if el_id in vdata:
                            st.annotations['rej' + fc] = True
                        else:
                            st.annotations['rej' + fc] = False
                    elif self.metadata_dir is not None:
                        # Annotate connector and brain aligned id
                        st.annotations['ca_id'] = \
                            block.annotations['elid_list_ca'].index(el_id) + 1

                    # Annotate if unit_id corresponds to sua or mua
                    if st.annotations['unit_id'] in self.get_sua_ids(el_id):
                        st.annotations['sua'] = True
                    else:
                        st.annotations['sua'] = False
                    if st.annotations['unit_id'] in self.get_mua_ids(el_id):
                        st.annotations['mua'] = True
                    else:
                        st.annotations['mua'] = False

        # Add annotations to units
        for unit in block.list_units:
            if unit.annotations['channel_id'] <= 100:
                el_id = unit.annotations['channel_id']
                unit.annotations['electrode_id'] = el_id
                if self.odML_avail:
                    # Annotate connector and brain aligned id
                    idx = block.annotations['elid_list'].index(el_id)
                    unit.annotations['ca_id'] = \
                        block.annotations['elid_list_ca'][idx]
                    st.annotations['ba_id'] = \
                        block.annotations['elid_list_ba'][idx]

                    # Annotate if electrode should be rejected
                    ff = lambda x: x.name == 'RejElectrodes' and \
                        'RejectionsLFP' in x.get_path()
                    pobjs = [p for p in
                             self.odML_doc.iterproperties(filter_func=ff)]
                    vdata = [v.data for v in pobjs[0].values]
                    fc = pobjs[0].parent.name
                    if el_id in vdata:
                        unit.annotations['rej' + fc] = True
                    else:
                        unit.annotations['rej' + fc] = False

                # Annotate if unit_id corresponds to sua or mua
                if unit.annotations['unit_id'] in self.get_sua_ids(el_id):
                    unit.annotations['sua'] = True
                else:
                    unit.annotations['sua'] = False
                if unit.annotations['unit_id'] in self.get_mua_ids(el_id):
                    unit.annotations['mua'] = True
                else:
                    unit.annotations['mua'] = False

        return block

    def add_trials_around(self, block, name, aligned_trigger, tpre, tpost,
                          aligned_offset=0 * pq.ms, performance_codes=[127]):
        """
        Adds trials of equal length as epocharray to the given neo block.

        Args:
            block (neo.block):
                neo Block object contatining the data.
            name (str):
                Name of the defined EpochArray
            aligned_trigger (int or str):
                Index or shorthand string of the event that should be used as
                trigger. See ReachGraspIO.trial_events for a list of event
                shorthands and ReachGraspIO.trial_events_codes for mapping them
                to corresponding indices, e.g. 2 for "TS".
            tpre (quantities.Quantity):
                Time that should be considered as trial before each trigger
                event.
            tpost (quantities.Ouantity):
                Time that should be considered as trial after each trigger
                event.
            align_offset (quantities.Ouantity):
                Time offset for time reset to aligned_trigger:
                    aligned_offset < 0: time reset point before aligned_trigger
                    aligned_offset > 0: time reset point after aligned_trigger
                    aligned_offset = 0: time reset point equals aligned_trigger
                Default: 0

        Returns:
            --

        Attributes and Annotations of EpochArray:
            EpochArray
                Per Segment, trials are stored in an EpochArray.
                Attributes:
                    name (str):
                        name of the defined EpochArray
                    file_origin (str):
                        session name
                    dtype:
                        depends on the units given for tpre and tpost
                    labels (list of str):
                        list of strings of the form "Trial ID] for each trial
                Annotations:
                    type (str):
                        equals 'trial' to identify these events as trial events
                    definition (str):
                        string defining trigger event, tpre and tpost
                    aligned_trigger (str):
                        shorthand for used trial_event
                    alignment (quantities.Quantity):
                        all time offsets for time reset to aligned_trigger
                    trial_id (list):
                        trial IDs
                    trial_type (list):
                        trial type codes
                    rej_mask_LFC (list):
                        mask for trial rejection for low frequency components
                        of the LFP (good trials == True)
                    rej_mask_IFC (list):
                        mask for trial rejection for intermediate frequency
                        components of the LFP (good trials == True)
                    rej_mask_HFC (list):
                        mask for trial rejection for high frequency components
                        of the LFP (good trials == True)
                    rej_mask_Spikes (list):
                        mask for trial rejection for spike data
        """
        # transform aligned_trigger from string to event index for trial data
        if type(aligned_trigger) == str:
            aligned_trigger = self.trial_data_idx[aligned_trigger]

        # transform performance codes to integer code if given as string
        if type(performance_codes) is not list:
            performance_codes = [performance_codes]
        for i, pc in enumerate(performance_codes):
            if type(pc) is str:
                performance_codes[i] = self.performance_code_idx[pc]

        # check if performance codes are compatible with given aligned trigger
        not_det_ev = []
        for pc in performance_codes:
            pc_binary = np.binary_repr(pc)
            for i, ev in enumerate([5, 6, 7, 8, 9, 11, 12]):
                if not int(pc_binary[i]):
                    not_det_ev.append(ev)
        if aligned_trigger in not_det_ev:
            raise IOError("Unable to generate trials with given performance" +
                          " code with the given aligned_trigger.")

        # Find trials in data
        aligned_trigger_times = \
            pq.Quantity(self.trial_data[:, aligned_trigger], self.nev_unit)

        # create mask for trials to consider according to performance_codes
        trial_ids = self.get_trial_ids()
        trial_ids_depr = self.get_trial_ids(trial_id_type='trial_id_depr')
        mask = (np.ones_like(trial_ids) == 1)
        if performance_codes:
            trial_ids = self.get_trial_ids(performance_codes=performance_codes)
            trial_ids_depr = \
                self.get_trial_ids(performance_codes=performance_codes,
                                   trial_id_type='trial_id_depr')
            mask = np.array([i in trial_ids for i in self.get_trial_ids()])

        # Get start times of all trials
        start_times = aligned_trigger_times[mask] - tpre
        # Get stop times of all trials
        stop_times = aligned_trigger_times[mask] + tpost
        # Get aligned_offsets of all trials
        aligned_offset_times = aligned_trigger_times[mask] + aligned_offset

        # Get trial types
        trial_types = np.ravel(
            self.trial_data[:, self.trial_data_idx['trial_type']])[mask]
        # Get performance code
        trial_pcs = np.ravel(
            self.trial_data[:, self.trial_data_idx['performance_code']])[mask]

        # Add trials as epocharray
        for seg in block.segments:
            # Find only those trials which are in the current segment
            idx = np.nonzero(
                np.logical_and(start_times >= seg.annotations['t_start'],
                               stop_times < seg.annotations['t_stop']))[0]

            labels = [name + " " + str(i) for i in np.ravel(trial_ids)[idx]]
            aligned_trig = self.trial_events_str[aligned_trigger]
            definition = "Trials triggered from " + aligned_trig + "-" + \
                         str(tpre) + " to " + aligned_trig + "+" + str(tpost)

            ep = neo.EpochArray(times=start_times[idx],
                                durations=stop_times[idx] - start_times[idx],
                                labels=labels,
                                name=name,
                                file_origin=self.associated_fileset,
                                type="trial",
                                definition=definition,
                                aligned_trigger=aligned_trig,
                                alignment=aligned_offset_times[idx],
                                trial_id=np.ravel(trial_ids)[idx],
                                trial_id_depr=np.ravel(trial_ids_depr)[idx],
                                trial_type=trial_types[idx],
                                trial_performance=trial_pcs[idx])

            # Add odML meta data info for each trial
            if self.odML_avail:
                trial_ids = self.get_trial_ids(performance_codes=127)
                # Annotate if electrode should be rejected
                ff = lambda x: x.name == 'RejTrials' and \
                    'RejectionsLFP' in x.get_path()
                pobjs = [p for p in
                         self.odML_doc.iterproperties(filter_func=ff)]
                rejtrials = [v.data for v in pobjs[0].values]
                fc = pobjs[0].parent.name
                ep.annotations['rej_mask_' + fc] = \
                    np.logical_not(np.in1d(trial_ids, rejtrials))

            seg.epocharrays.append(ep)
            seg.create_relationship()

    def add_trials_between(self, block, name, start_trigger, stop_trigger,
                           start_offset=0 * pq.ms, stop_offset=0 * pq.ms,
                           aligned_trigger=None, aligned_offset=0 * pq.ms,
                           performance_codes=[127]):
        """
        Adds trials of equal length as epocharray to the given neo block.

        Args:
            block (neo.block):
                neo Block object contatining the data.
            name (str):
                Name of the defined EpochArray
            start_trigger (int or str):
                Index or shorthand string of the event that should be used as
                start_trigger. See ReachGraspIO.trial_events for a list of event
                shorthands and ReachGraspIO.trial_events_codes for mapping them
                to corresponding indices, e.g. 2 for "TS".
            stop_trigger (int or str):
                Index or shorthand string of the event that should be used as
                stop_trigger. See ReachGraspIO.trial_events for a list of event
                shorthands and ReachGraspIO.trial_events_codes for mapping them
                to corresponding indices, e.g. 2 for "TS".
            start_offset (quantities.Ouantity):
                Time offset for the start_trigger:
                    start_offset < 0: time reset point before start_trigger
                    start_offset > 0: time reset point after start_trigger
                    start_offset = 0: time reset point equals start_trigger
                Default: 0
            stop_offset (quantities.Ouantity):
                Time offset for the stop_trigger:
                    stop_offset < 0: time reset point before stop_trigger
                    stop_offset > 0: time reset point after stop_trigger
                    stop_offset = 0: time reset point equals stop_trigger
                Default: 0
            aligned_trigger (int or str):
                Index or shorthand string of the event that should be used as
                aligned_trigger for the time reset. See ReachGraspIO.trial_events
                for a list of event shorthands and ReachGraspIO.trial_events_codes
                for mapping them to corresponding indices, e.g. 2 for "TS". If
                None is specified the aligned_trigger equals the start_trigger.
                Default: None
            align_offset (quantities.Ouantity):
                Time offset for time reset to aligned_trigger:
                    aligned_offset < 0: time reset point before aligned_trigger
                    aligned_offset > 0: time reset point after aligned_trigger
                    aligned_offset = 0: time reset point equals aligned_trigger
                Default: 0

        Returns:
            --

        Attributes and Annotations of EpochArray:
            EpochArray
                Per Segment, trials are stored in an EpochArray.
                Attributes:
                    name (str):
                        name of the defined EpochArray
                    file_origin (str):
                        session name
                    dtype:
                        depends on the units given for tpre and tpost
                    labels (list of str):
                        list of strings of the form "Trial ID] for each trial
                Annotations:
                    type (str):
                        equals 'trial' to identify these events as trial events
                    definition (str):
                        string defining trigger event, tpre and tpost
                    aligned_trigger (str):
                        shorthand for used trial_event
                    alignment (quantities.Quantity):
                        all time offsets for time reset to aligned_trigger
                    trial_id (list):
                        trial IDs
                    trial_type (list):
                        trial type codes
                    rej_mask_Low (list):
                        mask for trial rejection for low LFP frequencies
                    rej_mask_Mid (list):
                        mask for trial rejection for mid LFP frequencies
                    rej_mask_High (list):
                        mask for trial rejection for high LFP frequencies
                    rej_mask_Spikes (list):
                        mask for trial rejection for spike data
        """
        # get event index for start_trigger
        if type(start_trigger) is str:
            start_trigger = self.trial_data_idx[start_trigger]
        # get event index for stop_trigger
        if type(stop_trigger) is str:
            stop_trigger = self.trial_data_idx[stop_trigger]
        # get event index for aligned_trigger
        if type(aligned_trigger) is str:
            aligned_trigger = self.trial_data_idx[aligned_trigger]
        if aligned_trigger is None:
            aligned_trigger = start_trigger

        # transform performance codes to integer code if given as string
        if type(performance_codes) is not list:
            performance_codes = [performance_codes]
        for i, pc in enumerate(performance_codes):
            if type(pc) is str:
                performance_codes[i] = self.performance_code_idx[pc]

        # check if performance codes are compatible with given aligned trigger
        not_det_ev = []
        for pc in performance_codes:
            pc_binary = np.binary_repr(pc)
            for i, ev in enumerate([5, 6, 7, 8, 9, 11, 12]):
                if not int(pc_binary[i]):
                    not_det_ev.append(ev)
        if start_trigger in not_det_ev or stop_trigger in not_det_ev or \
           aligned_trigger in not_det_ev:
            raise IOError("Unable to generate trials with given performance" +
                          " code with the given aligned_trigger.")

        # create mask for trials to consider according to performance_codes
        trial_ids = self.get_trial_ids()
        trial_ids_depr = self.get_trial_ids(trial_id_type='trial_id_depr')
        mask = (np.ones_like(trial_ids) == 1)
        if performance_codes:
            trial_ids = self.get_trial_ids(performance_codes=performance_codes)
            trial_ids_depr = \
                self.get_trial_ids(performance_codes=performance_codes,
                                   trial_id_type='trial_id_depr')
            mask = np.array([i in trial_ids for i in self.get_trial_ids()])
        if len(mask) > 0:
            # Find trigger_times for all trials in trial_data
            start_trigger_times = \
                pq.Quantity(self.trial_data[:, start_trigger][mask],
                            self.nev_unit)
            stop_trigger_times = \
                pq.Quantity(self.trial_data[:, stop_trigger][mask],
                            self.nev_unit)
            aligned_trigger_times = \
                pq.Quantity(self.trial_data[:, aligned_trigger][mask],
                            self.nev_unit)
        else:
            raise ValueError("There are no trials with the given performance " +
                             "code detected which contain also the given " +
                             "start and stop trigger.")

        if len(start_trigger_times) != len(stop_trigger_times):
            raise ValueError("Length of start_times must equal length of " +
                             "stop_times")

        if False in [start_trigger_times[i] < stop_trigger_times[i] for i in
                     range(len(start_trigger_times))]:
            raise ValueError("start_times must be smaller than stop_times")

        # Get trial types
        trial_types = np.ravel(
            self.trial_data[:, self.trial_data_idx['trial_type']])[mask]
        # Get performance code
        trial_pcs = np.ravel(
            self.trial_data[:, self.trial_data_idx['performance_code']])[mask]

        # Add trials as epocharray
        for seg in block.segments:
            # Find only those trials which are in the current segment
            idx = np.nonzero(np.logical_and(
                start_trigger_times >= seg.annotations['t_start'],
                stop_trigger_times < seg.annotations['t_stop']))[0]

            durations = stop_trigger_times[idx] - start_trigger_times[idx]
            labels = [name + " " + str(i) for i in np.ravel(trial_ids)[idx]]
            start_trig = self.trial_events_str[start_trigger]
            stop_trig = self.trial_events_str[stop_trigger]
            aligned_trig = self.trial_events_str[aligned_trigger]
            definition = "Trials triggered from " + start_trig + "+" + \
                         str(start_offset) + " to " + stop_trig + "+" + \
                         str(stop_offset)

            ep = neo.EpochArray(times=start_trigger_times[idx],
                                durations=durations,
                                labels=labels,
                                name=name,
                                file_origin=self.associated_fileset,
                                type="trial",
                                defintion=definition,
                                start_trigger=start_trig,
                                stop_trigger=stop_trig,
                                start_offset=start_offset,
                                stop_offset=stop_offset,
                                aligned_trigger=aligned_trig,
                                alignment=aligned_offset,
                                trial_id=np.ravel(trial_ids)[idx],
                                trial_id_depr=np.ravel(trial_ids_depr)[idx],
                                trial_type=trial_types[idx],
                                trial_performance=trial_pcs[idx])

            # Add odML meta data info for each trial
            if self.odML_avail:
                trial_ids = np.ravel(self.trial_data[:, 0])
                for fr in ['Low', 'Mid', 'High']:
                    ff = lambda x: x.name == 'unique' and \
                        'Trials/Rejections/' + fr in x.parent.get_path()
                    vobjs = [p.values
                        for p in self.odML_doc.iterproperties(filter_func=ff)]
                    rejtrials = [v.data for v in vobjs[0]]
                    if type(rejtrials) is int and rejtrials < 0:
                        pass
                    elif type(rejtrials) is int and rejtrials > 0:
                        ep.annotations['rej_mask_' + fr] = \
                            np.logical_not(np.in1d(trial_ids, [rejtrials]))
                    else:
                        ep.annotations['rej_mask_' + fr] = \
                            np.logical_not(np.in1d(trial_ids, rejtrials))
            seg.epocharrays.append(ep)


def keep_synchrofacts(block, n=2, dt=0, dt2=1, unit_type=['sua']):
    """
    Keeps only spike artefacts in spiketrains of given unit type in given block.
    Accepted unit types: 'sua','mua','idx' (whereas x is the id number requested)
    unit_type = [] selects all available spike trains
    """
    del_spikes(block, n=n, dt=dt, dt2=dt2, invert=False, unit_type=unit_type)


def remove_synchrofacts(block, n=2, dt=0, dt2=1, unit_type=['sua']):
    """
    Removes all spike artefacts in spiketrains of block.
    Accepted unit types: 'sua','mua','idx' (whereas x is the id number requested)
    unit_type = [] selects all available spike trains
    """
    del_spikes(block, n=n, dt=dt, dt2=dt2, invert=True, unit_type=unit_type)


# _____________ Helper Functions for synchrofact handling ________________
def del_spikes(block, n, dt, dt2, invert=True, unit_type=['sua']):
    '''
    Given block with spike trains, delete all spikes engaged
    in synchronous events of size *n* or higher. If specified, delete
    spikes close to such syncrhonous events as well.

    *Args*
    ------
    block [list]:
        a block containing neo spike trains

    n [int]:
        minimum number of coincident spikes to report synchrony

    dt [int. Default: 0]:
        size of timestep lag for synchrony. Spikes closer than *dt* are
        considered synchronous. Groups of *n* or more synchronous spikes
        are deleted from the spike trains.

    dt2 [int. Default: 1]:
        maximum distance of timesteps for two spikes to be "close". Spikes
        "close" to synchronous spikes are eliminated as well.

    invert [bool. Default: True]:
        selects synchronous events to be deleted (Default:True). False keeps only
        synchrofacts and definded neigboring spikes in spiketrains

    unit_type [list of strings. Default ['sua']]:
        selects only spiketrain of certain units / channels for synchrofact extraction.
        unit_type = [] considers all provided spiketrains
        Accepted unit types: 'sua','mua','idx' (whereas x is the id number requested)
        Warning: for invert=True non-selected spiketrains are emptied, for invert=False
                no change occurs in non-selected spiketrains
    '''

    # data check
    if len(block.segments[0].spiketrains) == 0:
        warnings.warn(
            'Attempted synchronous event extraction with empty block')
        return block

    dt = dt * block.segments[0].spiketrains[0].units
    dt2 = dt2 * block.segments[0].spiketrains[0].units

    # extracting spiketrains which should be used for synchrofact extraction based on given unit type
    # possible improvement by using masks for different conditions and adding
    # them up
    neo_spiketrains, index = [], []
    for i in range(len(block.segments[0].spiketrains)):
        take_it = False
        for utype in unit_type:
            if utype[:2] == 'id' and block.segments[0].spiketrains[i].annotations['unit_id'] == int(utype[2:]):
                take_it = True
            elif (utype == 'sua' or utype == 'mua') and utype in block.segments[0].spiketrains[i].annotations and block.segments[0].spiketrains[i].annotations[utype]:
                take_it = True
        if take_it:
            neo_spiketrains.append(
                copy.deepcopy(block.segments[0].spiketrains[i]))
            index.append(i)

    # considering all spiketrains for unit_type == []
    if unit_type == []:
        neo_spiketrains = copy.deepcopy(block.segments[0].spiketrains)
        index = range(len(block.segments[0].spiketrains))

    # if no spiketrains were selected
    if len(neo_spiketrains) == 0:
        warnings.warn(
            'No matching spike trains for given unit selection criteria %s found' % (unit_type))
        times, dt = np.array([]), 0 * pq.s
        time_dt2 = np.array([])
    else:
        # find times of synchrony of size >=n
        times, dt = detect_syn_spikes(neo_spiketrains, n=n, dt=dt)
        # hstack only works for input != []
        times_dt2 = np.hstack([np.arange(
            j - dt2.base, j + dt2.base + 1) for j in times.magnitude]) if times != [] else []

    j = 0  # index of pre selected sts
    # iterating over original spike trains
    for idx in range(len(block.segments[0].spiketrains)):
        # for j, idx in enumerate(index):  # and copy in it the original ones
        # devoided of the synchrony times
        if idx in index:

            annotations = block.segments[0].spiketrains[idx].annotations

            mask = np.in1d(neo_spiketrains[j], times_dt2)
            if invert:
                mask = np.invert(mask)

            block.segments[0].spiketrains[
                idx] = neo_spiketrains[j][mask.nonzero()]

            if neo_spiketrains[j].waveforms != None:
                block.segments[0].spiketrains[idx].waveforms = neo_spiketrains[
                    j].waveforms[mask.nonzero()]

            block.segments[0].spiketrains[idx].annotations = annotations
            j += 1
        # case st was not selected -> removal of spikes in case of 'keep
        # synchrofacts'
        else:
            if not invert:
                block.segments[0].spiketrains[idx] = neo.SpikeTrain(np.array([]) * pq.ms,
                                                        t_start=block.segments[
                                                            0].spiketrains[idx].t_start,
                    t_stop=block.segments[
                                                            0].spiketrains[idx].t_stop,
                    **block.segments[0].spiketrains[idx].annotations)


# self, x, t_stop=10, t_start=0):
def spiketrains2gdf(neo_spiketrains, ids=[]):
    """
    Converts a list of spike trains to gdf format.

    Gdf is a 2-column data structure containing neuron ids on the first
    column and spike times (sorted in increasing order) on the second
    column. Information about the time unit, not preserved in the float-
    like gdf, is returned as a second output

    *Args*
    ------
    sts [list]:
        a list of neo spike trains.

    ids [list. Default to []]:
        List of neuron IDs. Id[i] is the id associated to spike train
        sts[i]. If empty list provided (default), dds are assigned as
        integers from 0 to n_spiketrains-1.

    *Returns*
    ---------
    gdf [ndarray of floats with shape (n_spikes, 2)]:
        ndarray of unit ids (first column) and
    """

#        t_stop = max([train.t_stop for train in neo_spiketrains])
#        t_start = min([train.t_start for train in neo_spiketrains])

    # Find smallest time unit
    time_unit = neo_spiketrains[0].units
    for st in neo_spiketrains[1:]:
        if st.units < time_unit:
            time_unit = st.units

#        self.unit = time_unit
#        self.times = neo.SpikeTrain(self[:, 1] * pq.second, t_start=t_start, t_stop=t_stop)
#        self.dtypes = (int, neo.SpikeTrain)

    # By default assign integers 0,1,... as ids of sts[0],sts[1],...
    if len(ids) == 0:
        ids = range(len(neo_spiketrains))

    gdf = np.zeros((1, 2))
    # Rescale all spike trains to that time unit, and add to the gdf
    for st_idx, st in zip(ids, neo_spiketrains):
        to_be_added = np.array([[st_idx] * len(st),
            st.view(pq.Quantity).rescale(time_unit).magnitude]).T
        gdf = np.vstack([gdf, to_be_added])

    # Eliminate first row in gdf and sort the others by increasing spike
    # times
    gdf = gdf[1:]
    gdf = gdf[np.argsort(gdf[:, 1])]

    # Return gdf and time unit corresponding to second column
    return gdf, time_unit


def detect_syn_spikes(neo_spiketrains, n=2, dt=0 * pq.ms, ids=[]):
    """
    Given a list *sts* of spike trains, returns the times where *n* or more
    spike trains have exactly synchronous spikes, and the neuron ids
    firing at those times

    *Args*
    ------
    sts [list]:
        a list of neo spike trains

    n [int]:
        minimum number of coincident spikes to report synchrony

    dt [Quantity. Default to 0 ms]:
        size of time lag for synchrony. Starting from the very first spike,
        spike trains are binned at a binsize *dt* (bins half-open to the
        right), and spikes in the same bin are considered synchronous. If 0
        (default)

    ids [list. Default to []]:
        List of neuron IDs. Id[i] is the id associated to spike train
        sts[i]. If empty list provided (default), ids are assigned as
        integers from 0 to n_spiketrains-1.


    *Returns*
    ---------
    tracts [list of ndarrays]:
        a list of transactions, each as an array of neuron ids

    times [Quantity]
        a list of coincidence times. coinc[i] is the time of tracts[i]
    """
    # Convert list of spike trains to time-sorted gdf
    gdf, time_unit = spiketrains2gdf(neo_spiketrains, ids=ids)
    dt_dimless = dt.rescale(time_unit).magnitude  # make dt dimension-free
    # if dt_dimless is 0, set to half the minimum non-zero ISI
    if dt_dimless == 0:
        dt_dimless = np.diff(np.unique(gdf[:, 1])).min() / 2.

    # TODO: Clean up comments
    # tracts, times = [], []  # Initialize transactions and times to be
    # returned to empy lists
    idx_synch = []  # numpy.empty((0,2))
    # Set the initial time for the synchrony search to first available spike
    # time
    time = gdf[0, 1]
    # idx_start, idx_stop = 0, 0 # Starting from the first row in the gdf
    idx_start, idx_stop = 0, 0  # starting from the very first spike in the gdf
    while idx_stop <= gdf.shape[0] - 2:  # until end of gdf is reached,
        # Until the next spike falls in [time, time+dt] (symmetrical!)
        while time <= gdf[idx_stop + 1, 1] <= time + dt_dimless:
            idx_stop += 1  # Include that spike in the transaction
            if idx_stop >= gdf.shape[0] - 1:
                break  # And stop if end of gdf reached
        # If at least n spikes fall between idx_start and idx_stop
        if idx_stop >= idx_start + n - 1:
            idx_synch.extend(range(idx_start, idx_stop + 1))
        idx_start += 1  # Set new idx_start to the next spike
        idx_stop = idx_start  # and idx_stop to idx_start
        time = gdf[idx_stop, 1]  # and set the new corresponding spike time.
    # check if last entries also fulfill synchrony condition
    # If at least n spikes fall between idx_start and idx_stop
    if idx_stop >= idx_start + n - 1:
        idx_synch.extend(range(idx_start, idx_stop + 1))

    idx_synch = np.array(np.unique(idx_synch), dtype=int)

    # Return transactions of >=n synchronous spikes, and the times of these
    # transactions (first spike time in each transaction)
    # #return tracts, times*time_unit, dt_dimless*time_unit
    return gdf[idx_synch][:, 1] * time_unit, dt_dimless * time_unit
