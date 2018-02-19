"""
Various tools to select SpikeTrain and AnalogSignal objects from a neo structure
"""

import numpy as np
import neo.core
import quantities as pq
from warnings import warn


def get_trials_AnalogSignal(bl, trial_prop, signal_prop, reset_times=True, cut_in_range=False):
    """
    This function returns a tuple with lists of AnalogSignal Objects, each corresponding to one trial,
    and corresponding EpochArray Object.

    Parameter:
    ----------
    bl : neo.Block
        The neo Block object to extract data from.
    trial_prop : dictionary
        A dictionary that contains selection criteria for trials. Each
        key of the dictionary is matched to an annotation of the
        EpochArray that contains a list of values for a certain property
        of that trial (eg., the trial ID). The value associated with
        the key contains the set of allowed values, i.e., which trials
        are to be considered.
        Note: If an epoch does not have a specific annotation, it is rejected.
    signal_prop : dictionary
        Same as trial_prop, but for annotations of the AnalogSignal object.
    reset_times : bool
        If True (t_start, t_stop) of analog signals are set to (0 *pq.s, epoch duration)
        If False (t_start, t_stop) are kept at original values.
        If reset_times is a string the annotation of the epoch with this string as
        keyword is looked up and the list of times is used as time points zero. (deprecated)
        Default: True
    cut_in_range : bool
        False: If the trial duration is longer than or exceeds the signal, the AnalogSignal
        Object will not be considered and appended to the list. (default)
        True: Even if the trial duration is longer than the signal and the staring point of the
        signal lies somewhere in the given range than the AnalogSignal Object will be considered and appended to
        the list.

    Returns:
    ---------
    list of tuples:
        A tuple of AnalogSignal objects (list) per trial and corresponding EpochArray.
        [(EpochArray, [AnalogSignalObject, ...])]

    Raises:
    --------
        TypeError: An error occurs if a wrong type is given for the parameter reset_time

    See also:
        get_trials_SpikeTrain

    """
    # List of AnalogSignal objects
    asig_lst = []
    # List of returned AnalogSignal objects in combination with EpochArrays, stored as tuple
    r = []

    # Step through all segments of the Block separately
    for seg in bl.segments:
        epochs = None
        for ep in seg.epocharrays:
            epochs = ep
            # TODO: hardcoding this is bad
            if ep.annotations['type'] == 'trial':
                # By default take all trials
                take_trial = np.array([True for _ in xrange(len(ep.durations))])

                # Now remove trials based on user _input
                for k in trial_prop.keys():
                    take_trial = take_trial & np.array([_ in trial_prop[k] for _ in ep.annotations[k]])

                # Step through all AnalogSignal objects
                for fp in seg.analogsignals:
                    take_analogsignal = True
                    for k in signal_prop.keys():
                        take_analogsignal = take_analogsignal and (fp.annotations[k] in signal_prop[k])

                    if take_analogsignal:
                        # Step through all trials
                        for i, t in enumerate(take_trial):
                            if t:
                                if _is_in_range(ep.times[i], ep.durations[i], fp.t_start, fp.t_stop,
                                                cut_in_range):
                                    # Range between edges
                                    ind1 = int(((ep.times[i] - fp.t_start) / (fp.t_stop - fp.t_start)).rescale('dimensionless') * len(fp))
                                    ind2 = int((
                                        ind1 + ep.durations[i].rescale(fp.t_start.units) / fp.sampling_period).base)
                                    # Don't trespass edge
                                    if ind2 > fp.t_start.magnitude + len(fp):
                                        ind2 = fp.t_start.magnitude + len(fp)

                                    # setting zero time to beginning of analogsignal
                                    if type(reset_times) == bool and reset_times is True:
                                        (t_start, t_stop) = (0 * ep.durations[i].units, ep.durations[i])

                                    # kept time stamps of original analogsignal
                                    elif type(reset_times == bool) and reset_times is False:
                                        (t_start, t_stop) = (ep.times[i], ep.times[i] + ep.durations[i])

                                        # setting zero time to provided event
                                        # elif type(reset_times) == str:
                                        #     if ep.annotations.has_key(reset_times):
                                        #         zero_times = ep.annotations[reset_times]
                                        #         if len(zero_times) == len(seg.epocharrays):
                                        #             # i corresponds to the number of current epoch
                                        #             (t_start, t_stop) = (ep.times[i] - zero_times[i], ep.times[i] - zero_times[i] + ep.durations[i])
                                        #
                                        #     else:
                                        #         raise ValueError("Provided keyword %s is not present in epoch annotations" % (reset_times))

                                    else:
                                        raise TypeError("Provided data type %s for reset times is not supported" % (
                                            type(reset_times)))
                                    asig_lst.append(neo.core.AnalogSignal(fp[ind1:ind2],
                                                                          t_stop=t_stop,
                                                                          sampling_rate=fp.sampling_rate,
                                                                          units=fp.units,
                                                                          t_start=t_start,
                                                                          trial_id=ep.annotations['trial_id'][i],
                                                                          **fp.annotations))
        if asig_lst and epochs:
            r.append((epochs, asig_lst))
    return r


def _is_in_range(ep_times, ep_durations, fp_start, fp_stop, cut_half):
    """
    Helper Function for get_trials, checks if cutting in (some) range is possible

    Parameter:
    ----------
    ep_times : int
        Starting time of neo.Epoch object.
    ep_durations : int
        Duration of neo.Epoch object.
    fp_start : int
        Starting time of signal.
    fp_stop : int
        Stopping time of signal
    cut_half : bool
        False: If the trial duration is longer than or exceeds the signal, the AnalogSignal
        Object will not be considered and appended to the list. (default)
        True: Even if the trial duration is longer than the signal and the staring point of the
        signal lies somewhere in the given range than the AnalogSignal Object will be considered and appended to
        the list.

    Returns:
    --------
    bool:
        True, if the cutting in given range is possible
        False, otherwise.
    """
    if cut_half:
        if fp_start <= ep_times <= fp_stop or fp_start <= ep_times + ep_durations <= fp_stop:
            return True
    else:
        if ep_times >= fp_start and ep_times + ep_durations <= fp_stop:
            return True
    return False


def get_trials_SpikeTrain(bl, trial_prop, signal_prop, reset_times=True):
    """
    This function returns a tuple with list of SpikeTrain Objects, each corresponding to one trial,
    and corresponding EpochArray Object.

    Parameter:
    ----------
    bl : neo.Block
        The neo Block object to extract data from.
    trial_prop : dictionary
        A dictionary that contains selection criteria for trials. Each
        key of the dictionary is matched to an annotation of the
        EpochArray that contains a list of values for a certain property
        of that trial (eg., the trial ID). The value associated with
        the key contains the set of allowed values, i.e., which trials
        are to be considered.
        Note: If an epoch does not have a specific annotation, it is rejected.
    signal_prop : dictionary
        Same as trial_prop, but for annotations of the SpikeTrain object.
    reset_times : bool
        If True (t_start, t_stop) of spiketrains are set to (0 *pq.s, epoch duration)
        If False (t_start, t_stop) are kept at original values.
        If reset_times is a string the annotation of the epoch with this string as
        keyword is looked up and the list of times is used as time points zero. (deprecated)
        Default: True

    Returns:
    --------
    list of tuples:
        A list of tuple of SpikeTrain objects (list) per trial and corresponding EpochArray.
        [(EpochArray, [SpikeTrainObject, ...])]

    Raises:
    ---------
    TypeError:
        An error occurs if a wrong type is given for the parameter reset_time

    See also:
        get_trials_AnalogSignal

        Examples:
            >>> get_trials_SpikeTrain(block,{trial_id=[1,2,3]})
            Returns a list of SpikeTrain objects corresponding to trials 1-3


            >>> get_trials_SpikeTrain(block,{trial_id=range(50),trial_type=[15,16]},...
                                      signal_prop={channel_id=15})

            Returns a list of SpikeTrain objects corresponding to trial type 15
            and 16, and where the trial ID is between 0 and 49, and recorded on
            channel 15.
    """
    # List of SpikeTrain objects
    spk_lst = []
    # List of returned SpikeTrain objects in combination with EpochArrays, stored as tuple
    r = []

    # Step through all segments of the Block separately
    for seg in bl.segments:
        epochs = None
        for ep in seg.epocharrays:
            epochs = ep
            # TODO: hardcoding this is bad
            if ep.annotations['type'] == 'trial':
                # By default take all trials
                take_trial = np.array([True for _ in xrange(len(ep.durations))])

                # Now remove trials based on user _input
                for k in trial_prop.keys():
                    take_trial = take_trial & np.array([_ in trial_prop[k] for _ in ep.annotations[k]])
                # Step through all SpikeTrain objects
                for st in seg.spiketrains:
                    take_spiketrain = True
                    for k in signal_prop.keys():
                        take_spiketrain = take_spiketrain and (st.annotations[k] in signal_prop[k])

                    if take_spiketrain:
                        # Step through all trials
                        # TODO: check if take_trial includes actually a trial
                        for i, t in enumerate(take_trial):
                            if t:
                                # setting zero time to beginning of spiketrain
                                if type(reset_times) == bool and reset_times is True:
                                    (t_start, t_stop) = (0 * st.units, ep.durations[i])
                                    spikingtimes = st.times[
                                        np.logical_and(
                                            st.times >= ep.times[i],
                                            st.times < ep.times[i] + \
                                            ep.durations[i])] - ep.times[i]

                                # kept time stamps of original spiketrain
                                elif type(reset_times == bool) and reset_times is False:
                                    (t_start, t_stop) = (ep.times[i], ep.times[i] + ep.durations[i])
                                    spikingtimes = st.times[
                                        np.logical_and(
                                            st.times >= ep.times[i],
                                            st.times < ep.times[i] + ep.durations[i])]

                                # setting zero time to provided event
                                # elif type(reset_times) == str:
                                #     if ep.annotations.has_key(reset_times):
                                #         zero_times = ep.annotations[reset_times]
                                #         if len(zero_times) == len(ep.times):
                                #             # i corresponds to the number of current epoch
                                #             (t_start, t_stop) = (ep.times[i] - zero_times[i], ep.times[i] - zero_times[i] + ep.durations[i])
                                #             spikingtimes = st.times[np.nonzero(np.logical_and(st.times >= ep.times[i], st.times < ep.times[i] + ep.durations[i]))] - zero_times[i]
                                #     else:
                                #         raise ValueError("Provided keyword %s is not present in epoch annotations" % (reset_times))
                                else:
                                    raise TypeError(
                                        "Provided data type %s for reset times is not supported" % (type(reset_times)))

                                if st.waveforms is not None:
                                    trial_waveforms = st.waveforms[np.nonzero(np.logical_and(st.times >= ep.times[i],
                                                                                             st.times < ep.times[i] +
                                                                                             ep.durations[i]))]

                                # [m for m in spikingtimes.magnitude] *st.units
                                trial_spiketrain = neo.core.SpikeTrain(
                                    pq.Quantity(spikingtimes.magnitude, units=spikingtimes.units),
                                    t_stop=t_stop,
                                    units=st.units,
                                    t_start=t_start,
                                    trial_id=ep.annotations['trial_id'][i],
                                    **st.annotations)
                                if st.waveforms is not None:
                                    trial_spiketrain.waveforms = trial_waveforms
                                spk_lst.append(trial_spiketrain)
        if epochs and spk_lst:
            r.append((epochs, spk_lst))
    return r


def get_sig_spiketrain(bl, signal_prop=None):
    """
    This function returns a list of SpikeTrain Objects, each corresponding to given id

    Parameter:
    ----------
    bl : neo.Block:
        The neo Block object to extract data from.
    signal_prop : dictionary
        A dictionary that contains the id numbers. For annotations of the SpikeTrain object.
        Each key of the dictionary is matched to an annotation of the id of the SpikeTrain object.
        Note: If no id or an empty dictionary is given all SpikeTrain Objects will be returned in a list.

    Returns:
    ---------
    list:
        A list of SpikeTrain objects.

    See also:
        get_sig_analogsignal

    Examples:
        >>> get_sig_spiketrain(block, id={1,2,3})
        Returns a list of SpikeTrain objects corresponding to id 1-3.
    """
    return _get_signal(bl, signal_prop, True)


def get_sig_analogsignal(bl, signal_prop=None):
    """
    This function returns a list of AnalogSignal Objects, each corresponding to given id.

    Parameter:
    ----------
    bl : neo.Block
        The neo Block object to extract data from.
    signal_prop : dictionary
        A dictionary that contains the id numbers. For annotations of the AnalogSignal object.
        Each key of the dictionary is matched to an annotation of the id of the AnalogSignal object.
         Note: If no id or an empty dictionary is given all AnalogSignal Objects will be returned in a list.

    Returns:
    ---------
    list:
        A list of AnalogSignal Objects.

    See also:
        get_sig_spiketrain

    Examples:
        >>> get_sig_analogsignal(block, id={1,2,3})
        Returns a list of AnalogSignal objects corresponding to id 1-3.
    """
    return _get_signal(bl, signal_prop, False)


def _get_signal(bl, signal_prop, is_spiketrain):
    # List of Signals
    sig_lst = []
    # SpikeTrain Routine
    if is_spiketrain:
        for seg in bl.segments:
            # If no or empty dictionary is given return all Spiketrains
            if not signal_prop or bool([b for b in signal_prop.values() if b == []]):
                sig_lst.append([spk for spk in seg.spiketrains])
            else:
                # Step through all SpikeTrain objects
                for spkt in seg.spiketrains:
                    take_spiketrain = True
                    for k in signal_prop.keys():
                        if k in spkt.annotations:
                            # Check if there is a key which resembles the given annotation
                            take_spiketrain = take_spiketrain and (spkt.annotations[k] in signal_prop[k])
                    if take_spiketrain:
                        sig_lst.append(spkt)
        return sig_lst
    # Analogsignal
    else:
        for seg in bl.segments:
            if not signal_prop or bool([b for b in signal_prop.values() if b == []]):
                sig_lst.append(asgl for asgl in seg.analogsignals)
            else:
                for asig in seg.analogsignals:
                    take_analogsignal = True
                    for k in signal_prop:
                        if k in asig.annotations:
                            take_analogsignal = take_analogsignal and (asig.annotations[k] in signal_prop[k])
                    if take_analogsignal:
                        sig_lst.append(asig)
        return sig_lst


def get_events(bl, event_prop=None):
    """
    This function returns a list of EventArray Objects, corresponding to given keywords.

    Parameter:
    ---------
    bl : neo.Block
        The neo Block object to extract data from.
    event_prop : dictionary
        A dictionary that contains the event properties.
        Each key of the dictionary is matched to an annotation of an object property.
        Note: If no properties or an empty dictionary is given, all EventArray Objects will be returned in a list.

    Returns:
    --------
        list:
            A list of EventArray Objects.
    """
    event_lst = []
    for seg in bl.segments:
        # empty or no dictionary
        if not event_prop or bool([b for b in event_prop.values() if b == []]):
            event_lst.append([e for e in seg.eventarrays])
        # dictionary is given
        else:
            for eva in seg.eventarrays:
                take_eva = True
                for k in event_prop.keys():
                    if k in eva.annotations:
                        take_eva = take_eva and (eva.annotations[k] in event_prop[k])
                        if take_eva:
                            event_lst.append(eva)
    return event_lst


def get_epochs(bl, epoch_prop=None):
    """
    This function returns a list of EpochArray Objects, corresponding to given keywords.

    Parameter:
    -------
    bl : neo.Block
        The neo Block object to extract data from.
    epoch_prop : dictionary
        A dictionary that contains the event properties.
        Each key of the dictionary is matched to an annotation of an object property.
        Note: If no properties or an empty dictionary is given, all EpochArray Objects will be returned in a list.

    Returns:
    -------
    list:
        A list of EpochArray Objects.
    """

    epoch_list = []
    for seg in bl.segments:
        # empty or no dictionary
        if not epoch_prop or bool([b for b in epoch_prop.values() if b == []]):
            epoch_list.append([e for e in seg.epocharrays])
        # dictionary is given
        else:
            for epa in seg.epocharrays:
                take_ep = True
                for k in epoch_prop.keys():
                    if k in epa.annotations:
                        take_ep = take_ep and (epa.annotations[k] in epoch_prop[k])
                        if take_ep:
                            epoch_list.append(epa)
    return epoch_list


def map_st(block, func, annotations=None, **kwargs):
    """
    The function 'func' with corresponding keywords 'kwargs' is applied to the SpikeTrain Objects of Block
    'block'.
    The function is applied in-place.

    Parameter:
    ----------
    bl : neo.Block
        The neo Block object to extract data (SpikeTrains) from.
    func:
        Function which will be applied on each SpikeTrain Object from Block
    annotations : dictionary
        User defined meta-data for provenance tracking.
        If **None** a String 'Functions' as key and the name of the function as value will be annotated.
        Raises a TypeError if a different data type than a dictionary is given.
    kwargs:
        To given function corresponding arguments as keywords

    Returns:
    ----------
        -

    Raises:
    TypeError, if a different data type than a dictionary is given.

    See also:
        map_as
    """
    for seg in block.segments:
        for idx, elem in enumerate(seg.spiketrains):
            seg.spiketrains[idx] = func(elem, **kwargs)
    if annotations is None:
        block.annotations["Functions: "] = func.__name__
    elif type(annotations) is dict and annotations:
        for key in annotations:
            block.annotations[key] = annotations[key]
    else:
        raise TypeError("Provided data type %s for annotations is not supported" % (type(annotations)))


def map_as(block, func, annotations, **kwargs):
    """
    The function 'func' with corresponding keywords 'kwargs' is applied to the AnalogSignal Objects
    of Block 'block'.
    The function is applied in-place.

    Parameter:
    ----------
    bl : neo.Block
        The neo Block object to extract data (AnalogSignals) from.
    func :
        Function which will be applied on each AnalogSignal Object from Block
    annotations : dictionary
        User defined meta-data for provenance tracking.
        If **None** a String 'Functions' as key and the name of the function as value will be annotated.
        Raises a TypeError if a different data type than a dictionary is given.
    kwargs:
        To given function corresponding arguments as keywords

    Returns:
    --------
        -

    Raises:
        TypeError, if a different data type than a dictionary is given.

    See also:
        map_st
    """
    for seg in block.segments:
        for idx, elem in enumerate(seg.analogsignals):
            seg.analogsignals[idx] = func(elem, **kwargs)
    if annotations is None:
        block.annotations["Functions: "] = func.__name__
    elif type(annotations) is dict and annotations:
        for key in annotations:
            block.annotations[key] = annotations[key]
    else:
        raise TypeError("Provided data type %s for annotations is not supported" % (type(annotations)))


def get_trials(bl, epoch_prop, as_prop=None, st_prop=None, reset_times=True,
        cut_in_range=False):
    """
    This function selects and cuts trial data from a neo Block object.

    The function returns a neo Block where each trial is represented as a
    neo.Segment in the Block, to which the cut SpikeTrain and/or
    AnalogSignal object are attached.

    Three dictionaries contain restrictions on which trials, and which data
    is included in the final data set. To this end, it is possible to specify
    accepted (valid) values of specific annotations on the source objects. Only
    those objects that match these restrictions are included in the final
    result.

    The resulting trials may either retain their original time stamps, or
    shifted to a common time axis.

    Parameters
    ----------
    bl : neo.Block
        The neo Block object to extract data from.
    epoch_prop : dictionary
        A dictionary that contains selection criteria for trials. Each
        key of the dictionary is matched to an annotation of the
        EpochArray that contains a list of values for a certain property
        of that trial (e.g., the trial ID). The value associated with
        the key contains the set of allowed values, i.e., which trials
        are to be considered.
        Note: If an epoch does not have a specific annotation, it is rejected.
    as_prop : dictionary
        Each key of as_prop is matched to an annotation of the AnalogSignal.
        The value associated to each key must be a list of allowed values for
        the annotation. Only those AnalogSignals are considered in the output,
        where (i) the corresponding annotation exists for the AnalogSignal and
        (ii) its value matches one of the values supplied v.
        If as_prop is None, no AnalogSignals are considered in the output
        Default: {}
    st_prop : dictionary
        Similar to as_prop, but for the SpikeTrain objects.
    reset_times : bool
        If True  the times stamps of all resulting object are set to fall in
            the range (0 *pq.s, epoch duration).
        If False, original time stamps are retained. Here the trials the
        Default is True
    cut_in_range : bool
        False: If the trial duration is longer than or exceeds the signal,
        the AnalogSignal will not be considered.
        True: If the trial duration is longer than a data object and tstart
            of the signal lies in the trial then the data object will be
            considered in the output.
        Default: False

    Returns:
    --------
    Block: neo.Block
        Per trial a neo.Segment with AnalogSignal or SpikeTrain Objects will
        be attached to a Block, if the properties of each Object are matching
        to the user given properties.
    """

    # TODO: Cannot deal with AnalogSignalArray and similar (although
    # soon redundant).
    if as_prop is None and st_prop is None:
        raise ValueError("No data specified to keep in the output.")
    block = neo.core.Block()
    for seg in bl.segments:
        for ep in seg.epocharrays:
            if ep.annotations['type'] == 'trial':
                # dictionary to hold all segments for valid trials
                segs = {}

                # by default take all trials
                take_trial = np.array(
                    [True for _ in xrange(len(ep.durations))])
                # now remove trials based on user _input
                for k in epoch_prop.keys():
                    take_trial = \
                        take_trial & \
                        np.array([
                            _ in epoch_prop[k] for _ in ep.annotations[k]])
                # step through all trials
                for i, t in enumerate(take_trial):
                    if t:
                        segs[i] = neo.core.Segment()

                # check if to include AnalogSignals
                if as_prop is not None:
                    # step through all AnalogSignal objects
                    for fp in seg.analogsignals:
                        take_analogsignal = True
                        for k in as_prop.keys():
                            take_analogsignal = \
                                take_analogsignal and \
                                (fp.annotations[k] in as_prop[k])
                        if take_analogsignal:
                            # step through all trials
                            for i, t in enumerate(take_trial):
                                if t:
                                    # check if in range
                                    if _is_in_range(
                                            ep.times[i], ep.durations[i],
                                            fp.t_start, fp.t_stop,
                                            cut_in_range):
                                        analog_sig = __cut_as(
                                            ep, fp, reset_times, i)

                                        # check if there is an AnalogSignal
                                        if np.size(analog_sig) > 0:
                                            segs[i].analogsignals.append(
                                                analog_sig)
                                        segs[i].create_relationship()
                # check if to include SpikeTrains
                if st_prop is not None:
                    # step through all SpikeTrain objects
                    for st in seg.spiketrains:
                        take_spiketrain = True
                        for k in st_prop.keys():
                            take_spiketrain = \
                                take_spiketrain and \
                                (st.annotations[k] in st_prop[k])
                        if take_spiketrain:
                            # Step through all trials
                            for i, t in enumerate(take_trial):
                                if t:
                                    if _is_in_range(ep.times[i],
                                                    ep.durations[i],
                                                    st.t_start, st.t_stop,
                                                    cut_in_range):
                                        st_sig = __cut_st(ep, st,
                                                          reset_times, i)
                                        # Check if there is a SpikeTrain
                                        if np.size(st_sig) > 0:
                                            segs[i].spiketrains.append(
                                                st_sig)
                                            segs[i].create_relationship()
                for s in sorted(segs):
                    block.segments.append(segs[s])
    block.create_relationship()
    return block


def __cut_as(ep, fp, reset_times, i):
    """
    Helper function to cut an AnalogSignal Object

    Parameter:
    ---------
    ep : neo.EpochArray
    fp : neo.AnalogSignal
    reset_times : bool
    i : int
        Actual index of the loop

    Returns:
    -------
    neo.AnalogSignal:
        A cut AnalogSignal Object.
    None:
        Only if the trial doesn't match inside the range and thus a AnalogSignal Object couldn't be constructed.
    """
    # Range between edges
    ind1 = int(((ep.times[i] - fp.t_start) / (fp.t_stop - fp.t_start)).rescale('dimensionless') * len(fp))
    ind2 = int((
        ind1 + ep.durations[i].rescale(fp.t_start.units) / fp.sampling_period).base)
    # Don't trespass edge
    if ind2 > fp.t_start.magnitude + len(fp):
        ind2 = fp.t_start.magnitude + len(fp)

    # setting zero time to beginning of analogsignal
    if type(reset_times) == bool and reset_times is True:
        (t_start, t_stop) = (0 * ep.durations[i].units, ep.durations[i])

    # keep time stamps of original spiketrain
    elif type(reset_times == bool) and reset_times is False:
        (t_start, t_stop) = (ep.times[i], ep.times[i] + ep.durations[i])
    else:
        raise TypeError("Provided data type %s for reset times is not supported" % (
            type(reset_times)))
    analog_sig = neo.core.AnalogSignal(fp[ind1:ind2],
                                       t_stop=t_stop,
                                       sampling_rate=fp.sampling_rate,
                                       units=fp.units,
                                       t_start=t_start,
                                       trial_id=ep.annotations['trial_id'][i],
                                       **fp.annotations)
    return analog_sig


def __cut_st(ep, st, reset_times, i):
    """
    Helper function to cut a SpikeTrain Object

    Parameter:
    ---------
    ep : neo.EpochArray
    fp : neo.AnalogSignal
    reset_times : bool
    i : int
        Actual index of the loop

    Returns:
    -------
    neo.SpikeTrain:
        A cut SpikeTrain Object.
    None:
        Only if the trial doesn't match inside the range and thus a SpikeTrain Object couldn't be constructed.
    """
    # setting zero time to beginning of spiketrain
    if type(reset_times) == bool and reset_times is True:
        (t_start, t_stop) = (0 * st.units, ep.durations[i])
        spikingtimes = st.times[np.nonzero(np.logical_and(st.times >= ep.times[i],
                                                          st.times < ep.times[i] +
                                                          ep.durations[i]))] - ep.times[i]
    # keep time stamps of original spiketrain
    elif type(reset_times == bool) and reset_times is False:
        (t_start, t_stop) = (ep.times[i], ep.times[i] + ep.durations[i])
        spikingtimes = st.times[np.nonzero(np.logical_and(st.times >= ep.times[i],
                                                          st.times < ep.times[i] +
                                                          ep.durations[i]))]

    else:
        raise TypeError(
            "Provided data type %s for reset times is not supported" % (type(reset_times)))

    if st.waveforms is not None:
        trial_waveforms = st.waveforms[np.nonzero(np.logical_and(st.times >= ep.times[i],
                                                                 st.times < ep.times[i] +
                                                                 ep.durations[i]))]

    if len(spikingtimes) < 1:
        warn('Empty SpikeTrain at trial %d' % i)
    trial_spiketrain = neo.core.SpikeTrain(
        pq.Quantity(spikingtimes.magnitude, units=spikingtimes.units),
        t_stop=t_stop,
        units=st.units,
        t_start=t_start,
        trial_id=ep.annotations['trial_id'][i],
        **st.annotations)
    if st.waveforms is not None:
        trial_spiketrain.waveforms = trial_waveforms
    return trial_spiketrain


