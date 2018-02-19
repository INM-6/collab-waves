# -*- coding: utf-8 -*-

import sys
import os
import math
import platform
import subprocess
import multiprocessing

import numpy as np
import quantities as pq
import scipy.cluster.hierarchy as hcluster
import scipy.cluster.vq as vqcluster
import scipy.signal as pysig

import neo

import pick


# =============================================================================
# Lists of datasets and global scan parameters
# =============================================================================


# List of monkey names and abbreviations of file names
monkey_names = {
    't': 'Tanya', 'a': 'Tanya2', 'l': 'Lilou',
    'n': 'Nikos', 'i': 'Nikos2', 's': 'Sana'}

# List of the 75 datasets and their condition code (1=grip-first, 2=force-
# first) considered in this project
selected_datasets = [
    ['l101210-001', 1],
    ['i140703-001', 1]]

# Short names for the used filters
selected_filter_names = [
    'beta_large',
    'beta_precise',
    'beta_high',
    'theta',
    'alpha',
    'gamma_low',
    'gamma_high']

# Parameters for each of the filters. If only one parameter set is given per
# filter, it is used for all monkeys. Otherwise, the dictionary associated to a
# filter lists the settings separately according to "monkey_names".
selected_filters = {
    selected_filter_names[0]: {
        'Tanya': {'lowcut': 13.0, 'highcut': 30.0, 'order': 4},
        'Tanya2': {'lowcut': 13.0, 'highcut': 30.0, 'order': 4},
        'Lilou': {'lowcut': 13.0, 'highcut': 30.0, 'order': 4},
        'Nikos': {'lowcut': 13.0, 'highcut': 30.0, 'order': 4},
        'Nikos2': {'lowcut': 13.0, 'highcut': 30.0, 'order': 4}},
    selected_filter_names[1]: {
        'Tanya': {'lowcut': 15.0, 'highcut': 21.0, 'order': 4},
        'Tanya2': {'lowcut': 15.0, 'highcut': 21.0, 'order': 4},
        'Lilou': {'lowcut': 19.0, 'highcut': 23.0, 'order': 4},
        'Nikos': {'lowcut': 15.0, 'highcut': 22.0, 'order': 4},
        'Nikos2': {'lowcut': 15.0, 'highcut': 22.0, 'order': 4}},
    selected_filter_names[2]: {'lowcut': 31.0, 'highcut': 48.0, 'order': 4},
    selected_filter_names[3]: {'lowcut': 4.0, 'highcut': 8.0, 'order': 2},
    selected_filter_names[4]: {'lowcut': 8.0, 'highcut': 13.0, 'order': 3},
    selected_filter_names[5]: {'lowcut': 52.0, 'highcut': 70.0, 'order': 5},
    selected_filter_names[6]: {'lowcut': 70.0, 'highcut': 100.0, 'order': 5}}

# Names of the 5 states
state_names = [
    'unclassified', 'planar', 'synchronized', 'random', 'circular',
    'radial']

state_names_abbrev = [
    'unclass.', 'planar', 'synchr.', 'random', 'circular',
    'radial']

# Condition names
cnd_names = {'1': 'grip-first', '2': 'force-first'}

# Degrees to turn array clockwise
array_turn = {
    'Lilou': 218. * 2. * np.pi / 360.,
    'Tanya': 320 * 2. * np.pi / 360.,
    'Nikos2': 239 * 2. * np.pi / 360.}

# =============================================================================
# Organizational functions
# =============================================================================


def fetch_data(sessionnames):
    """
    Determines whether the script is run on hambach. If not, then
    -- data of the sessions in sessionnames, and
    -- all metadata of all monkeys
    will be copied to the data cache on the laptop/blaustein in:
    ~/DatasetsCached/reachgrasp

    The subdirectory structure will be kept, such that data files reside in
    Data<Monkey> and metadata ind MetaData<Monkey>.

    Parameters
    ----------
    sessionnames : list of strings
        List of session names (without extensions, e.g., 'l101010-001') to copy
        to the cache, if on laptop.

    Note
    ----
    A call to this function will issue two calls to rsync to get data from the
    cluster. Calling this function often in sequence, or even worse, in a
    parallel fashion will bombard the server with ssh requests and lead to a
    security block to defend against what looks like a DDOS attack.
    """

    return

    # Nothing to do if we are on the server
    if platform.node()[:8] == 'hambach' or \
            platform.node()[:10] == 'blaustein':
        return

    # Directory for cached files on laptop
    data_dir = os.path.expanduser('~') + '/DatasetsCached/reachgrasp'

    # Initial command for copying data
    datacmd = ['rsync', '-a', '--progress', '--delete']

    # Include all Data Directories
    for k in monkey_names.values():
        datacmd.append('--include=Data' + k)

    # Exclude all mat files in general
    datacmd.append('--exclude=*.mat')

    # Include all sessions
    for k in sessionnames:
        # Append each file as include to command
        datacmd.append('--include=' + k + '*')

    # Exclude all other files
    datacmd.append('--exclude=*')

    # Append source to command
    datacmd.append(
        'denker@login.inm.kfa-juelich.de:' +
        '/datasets/marseille/congloue/data/DataGrasp/')

    # Append destination to command
    datacmd.append(data_dir + os.sep)

    # Initial command for copying metadata
    metadatacmd = ['rsync', '-a', '--progress', '--delete']

    # Copy all known metadata subdirectories
    for k in monkey_names.values():
        metadatacmd.append(
            # Replace with location of real odML files in time!
            'denker@login.inm.kfa-juelich.de:/datasets/' +
            # 'marseille/congloue/data/DataGrasp/MetaData' +
            'reachgrasp/metadata/MetaData' +
            k)

    metadatacmd.append(data_dir)

    # Now get the newest data and metadata files if not already available
    print(" ".join(datacmd))
    print(" ".join(metadatacmd))
    subprocess.call(datacmd)
    # subprocess.call(metadatacmd)


def get_data_dir(sessionname):
    """
    Determines whether the script is run on a laptop or on the cluster and
    returns the respective data directories accordingly.

    On the laptop or on blaustein, the cache is located at
    ~/DatasetsCached/reachgrasp

    Parameters
    ----------
    sessionname : string
        The sessionname to work on (without extensions, e.g., 'l101010').

    Returns
    -------
    data_dir : string
        Directory where to find the data (no trailing /).
    metadata_dir : string
        Directory where to find the metadata (no trailing /).
    """

    return ("datasets/","reachgrasp-spikewave/metadata/Meta" + subdir)

    # Determine subdirectory
    subdir = 'Data' + monkey_names[sessionname[0]]

    # Laptop or server
    if platform.node() != 'hambach':
        homedir = os.path.expanduser('~')
        # File will be on the external USB disk or stored locally on blaustein
        # as long as main data dir is not accessible
        data_dir = \
            homedir + '/DatasetsCached/reachgrasp/' + subdir
        metadata_dir = \
            homedir + '/DatasetsCached/reachgrasp/Meta' + subdir
    else:
        # Running on the server
        data_dir = \
            '/datasets/marseille/congloue/data/DataGrasp/' + subdir
        # Ignore odML files for now
#         metadata_dir = \
#             '/datasets/marseille/congloue/data/DataGrasp/Meta' + subdir
        metadata_dir = \
            '/datasets/reachgrasp/metadata/Meta' + subdir

    return (data_dir, metadata_dir)


def get_jobid():
    """
    For an embarrassingly trivial job, determines the job parameter
        - from PBS_ARRYAID if available, otherwise
        - from the first command-line argument, otherwise
        - a value of -1 to indicate that all jobs are to be launched

    Returns
    -------
    job_id: integer
        job_id to launch
    """

    # Get job parameter:
    # - PBS_ARRAYID
    # - or first input parameter
    # - or otherwise default: -1
    default_id = "-1"
    if len(sys.argv) > 1:
        cmdline_id = sys.argv[1]
        if cmdline_id.isdigit:
            default_id = cmdline_id
    job_id = int(os.getenv("PBS_ARRAYID", default_id))

    return job_id


def jobid2params(job_id, parameter_lists):
    """
    Convert a linear job ID to the indices of a given set of parameter lists.
    For example, for 3 parameter lists ['a','b','c'], [1,2,3,4], and [3.4,2.1],
    a job ID of 2 would result in indices (2,0,0), corresponding to parameters
    ['c',1,3.4]. A job ID of 7 would correspond to ['b',2,3.4]. The order of
    assigning job IDs to parameters if "Fortran" style, i.e., in our example,
    the first job IDs are mapped to
        (0,0,0)
        (1,0,0)
        (2,0,0)
        (0,1,0)
        (1,1,0)
        ...
        (2,3,0)
        (0,0,1)
        (1,0,1)
        ...

    Parameters
    ----------
        job_id: integer
            Linear job ID to convert.
        parameter_lists: list
            List of N lists of parameters to scan.

    Returns
    -------
        parameter_ids: list
            List of N indices into each of the parameter lists for the given
            job ID.

    Example
    -------
        >>> jobid2param(7, [['a', 'b', 'c'], [1, 2, 3, 4], [3.4, 2.1]])
        [1, 1, 0]
    """
    parameter_dims = []
    for p in parameter_lists:
        parameter_dims.append(len(p))
    return np.unravel_index(job_id, dims=parameter_dims, order='F')


def params2jobid(parameter_ids, parameter_lists):
    """
    Convert a list of indices into parameter lists into a linear job ID. The
    operation is the reverse of jobid2params. Please see the documentation of
    jobid2params for details.

    Parameters
    ----------
    parameter_ids: list
        List of N indices into each of the parameter lists.
    parameter_lists: list
        List of N lists of parameters to scan.

    Returns
    -------
    job_id: integer
        Linear job ID corresponding to the individual indices.

    Example
    -------
    >>> param2jobid([1, 1, 0], [['a', 'b', 'c'], [1, 2, 3, 4], [3.4, 2.1]])
    7
    """
    parameter_dims = []
    for p in parameter_lists:
        parameter_dims.append(len(p))
    return np.unravel_index(parameter_ids, dims=parameter_dims, order='F')


def create_file_filter_task(subsession_list, filter_list, fetch=True):
    """
    Creates a parameter list for tasks that require a combination of a file and
    a filter. The routine will also ensure that the data of all sub_sessions is
    available.

    Parameters
    ----------
    subsession_list: list
        A list of all subsessions (similar to selected_datasets).
    filter_list: list
        A list of all filters (similar to selected_filters).
    fetch: bool
        If True, fetch_data() is used to make sure all data files are
        available. Default: True

    Returns
    -------
    parameters: list
        List of lists of parameters for the job ID passed, or for all jobs if
        job_id was -1.
    """

    # Calculate number of jobs to run
    num_jobs = len(subsession_list) * len(filter_list)
    job_ids = range(num_jobs)

    # Fetch all data files (sequentially, to avoid DDOS attack on server!)
    parameters = []
    sessionnames = []
    for j in job_ids:
        param_indices = jobid2params(j, [subsession_list, filter_list])
        # Session to consider for the job
        selected_subsession_id = param_indices[0]
        selected_subsession = subsession_list[selected_subsession_id]

        # Filter range to consider for job
        selected_filter_id = param_indices[1]
        selected_filter = filter_list[selected_filter_id]

        print("Building job ID %i (%s,%s)" %
              (j, selected_subsession[0], selected_filter))

        # For each job, save its parameters
        parameters.append([selected_subsession, selected_filter])

        # Also save session names to fetch
        sessionnames.append(selected_subsession[0])

    # Pre-fetch data for this session if necessary
    if fetch:
        fetch_data(set(sessionnames))

    return parameters


def run_task(worker_function, job_id, parameter_list, spare_core=-1):
    """
    Launches either a single instance of the worker_function in batch queue
    mode, or uses multiprocessing to launch all jobs.

    Parameters
    ----------
    worker_function: object
        Function to be called as worker.
    job_id: integer or list
        Job IDs to execute. If an integer >=0, the corresponding job will be
        executed. If a list of integer job IDs, the multiprocessing package
        will be used to launch one job for each job ID. If -1, a job is
        launched for all parameters in parameter_list using the multiprocessing
        package.
    parameter_list: list
        List of lists of parameters to pass to the function for each job, i.e.,
        job ID 0 will be passed the contents of the list in parameter_list[0].
    spare_core: integer
        Number of cores to spare in multiprocessing. For example, if set to 2,
        only 6 cores on a 8 core machines will not be used. If -1, no
        multiprocessing is used at all. Default: -1

        Number of CPU cores to leave unoccupied in multiprocessing mode. At
        least one CPU must be available for running tasks, such that spare_core
        must be an integer smaller than the number of cores available on the
        machine. For example, if set to 2, only 6 cores on a 8 core machines 
        will not be used. If -1, no multiprocessing is used at all. 
        Default: -1
    """

    # If job_id is -1, use internally multiprocessing and process all sessions
    if job_id == -1:
        # We have as many job IDs as parameters in the parameter list
        job_id = range(len(parameter_list))

    if type(job_id) is list:
        # Set up multiprocessing with one CPU less than available
        if spare_core == -1:
            num_core = 1
        else:
            num_core = int(multiprocessing.cpu_count() - spare_core)

        if num_core < 1:
            raise ValueError(
                "The parameter spare_CPU must be smaller than the number of "
                "cores available on the current CPU")
        if num_core == 1:
            # No need for multiprocessing -- only one core (allows better
            # debugging)
            for j in job_id:
                worker_function(j, *parameter_list[j])
        else:
            # Use multiprocessing to parallelize
            pool = multiprocessing.Pool(processes=num_core)

            # Launch jobs
            results = []
            for j in job_id:
                results.append(pool.apply_async(
                    func=worker_function,
                    args=[j] + parameter_list[j]))
            pool.close()
            for r in results:
                r.get()
            print("Multiprocessing done.")
            pool.join()
            print("Multiprocessing joined.")
            pool.terminate()
    else:
        # Launch only a single job in batch queue mode
        worker_function(
            job_id, *parameter_list[job_id])


# =============================================================================
# Analytics
# =============================================================================


def norm_angle(p):
    """
    Maps an angle to the interval [-pi,pi).

    This function is useful to obtain the minimum angle between two anglular
    values a1 and a2 using:
        >>> da = norm_angle(a1-a2)
    In this case, mod(a1+da,2*pi)=mod(a2,2*pi) such that da<=pi, i.e., da is
    the smallest number to add to a1 such that the result is an angle in the
    same direction as a2.

    Parameters
    ----------
    angle: float
        Angle in rad.

    Returns
    -------
    norm_angle: float
        Normalized angle in (-pi,pi].
    """

    return(-np.mod(-p - np.pi, 2. * np.pi) + np.pi)

    # This is the solution for obtaining angles [-pi,pi):
    # return(np.mod(p + np.pi, 2. * np.pi) - np.pi)


def get_nn(electrode, nextneighbor_distance):
    """
    Returns a list of all next neighbors (independent of whether the electrode
    is broken or not). The electrode itself is NOT included.

    Parameters
    ----------
    electrode: int
        Linear electrode ID (1-100).
    nextneighbor_distance: int
        Number of electrodes in the horizontal, vertical, and diagonal
        direction that make up the neighborhood.

    Returns
    -------
    nn: list of int
        List of eletrode IDs in the neighborhood of the electrode.
    """

    nn = []

    el_row = int((electrode - 1) / 10) + 1
    el_col = np.mod(electrode - 1, 10) + 1

    for i in range(-nextneighbor_distance, nextneighbor_distance + 1):
        for j in range(-nextneighbor_distance, nextneighbor_distance + 1):
            if el_col + i >= 1 and el_col + i <= 10 and \
                    el_row + j >= 1 and el_row + j <= 10:
                newpos = electrode + i + j * 10
                if newpos != electrode:
                    nn.append(newpos)

    return sorted(nn)


def get_nn_nodiag(electrode, nextneighbor_distance):
    """
    Returns a list of all next neighbors (independent of whether the electrode
    is broken or not). The electrode itself is NOT included.

    Parameters
    ----------
    electrode: int
        Linear electrode ID (1-100).
    nextneighbor_distance: int
        Number of electrodes in the horizontal, vertical, and diagonal
        direction that make up the neighborhood.

    Returns
    -------
    nn: list of int
        List of eletrode IDs in the neighborhood of the electrode.
    """

    nn = []

    el_row = int((electrode - 1) / 10) + 1
    el_col = np.mod(electrode - 1, 10) + 1

    for i in range(-nextneighbor_distance, nextneighbor_distance + 1):
        if el_col + i >= 1 and el_col + i <= 10:
            newpos = electrode + i
            if newpos != electrode:
                nn.append(newpos)

        if el_row + i >= 1 and el_row + i <= 10:
            newpos = electrode + i * 10
            if newpos != electrode:
                nn.append(newpos)

    return sorted(nn)


def get_nn_row(electrode, nextneighbor_distance):
    """
    Returns a list of all next neighbors (independent of whether the electrode
    is broken or not). The electrode itself is NOT included.

    Parameters
    ----------
    electrode: int
        Linear electrode ID (1-100).
    nextneighbor_distance: int
        Number of electrodes in the horizontal, vertical, and diagonal
        direction that make up the neighborhood.

    Returns
    -------
    nn: list of int
        List of eletrode IDs in the neighborhood of the electrode.
    """

    nn = []

    el_col = np.mod(electrode - 1, 10) + 1

    for i in range(-nextneighbor_distance, nextneighbor_distance + 1):
        if el_col + i >= 1 and el_col + i <= 10:
            newpos = electrode + i
            if newpos != electrode:
                nn.append(newpos)

    return sorted(nn)


def get_nn_col(electrode, nextneighbor_distance):
    """
    Returns a list of all next neighbors (independent of whether the electrode
    is broken or not). The electrode itself is NOT included.

    Parameters
    ----------
    electrode: int
        Linear electrode ID (1-100).
    nextneighbor_distance: int
        Number of electrodes in the horizontal, vertical, and diagonal
        direction that make up the neighborhood.

    Returns
    -------
    nn: list of int
        List of eletrode IDs in the neighborhood of the electrode.
    """

    nn = []

    el_row = int((electrode - 1) / 10) + 1

    for i in range(-nextneighbor_distance, nextneighbor_distance + 1):
        if el_row + i >= 1 and el_row + i <= 10:
            newpos = electrode + i * 10
            if newpos != electrode:
                nn.append(newpos)

    return sorted(nn)


def get_elec_positions(electrode):
    """
    Returns a tuple corresponding to (X,Y) of the specified linear electrode
    ID.

    Parameters
    ----------
    electrode: int
        Linear electrode ID (1-100).

    Returns
    -------
    X: int
    Y: int
        X and Y position of the electrode.
    """

    return (np.mod(electrode - 1, 10) + 1, (electrode - 1) / 10 + 1)


def get_elec_distance(electrode1, electrode2):
    """
    Returns the distance between two electrodes.

    Parameters
    ----------
    electrode1: int
    electrode2: int
        Linear electrode IDs (1-100) of the two electrodes.

    Returns
    -------
    d: float
        Distance between the two electrodes in micrometers (um).
    """
    el1_pos = get_elec_positions(electrode1)
    el2_pos = get_elec_positions(electrode2)

    dist_row = (el1_pos[0] - el2_pos[0]) * 400
    dist_col = (el1_pos[1] - el2_pos[1]) * 400
    d = np.sqrt(dist_row ** 2 + dist_col ** 2)

    return d


def get_elec_angle(electrode1, electrode2):
    """
    Returns the angle between two electrodes.

    Parameters
    ----------
    electrode1: int
        Linear electrode IDs (1-100) of the first electrode from which the
        angle is measured.
    electrode2: int
        Linear electrode IDs (1-100) of the electrode to which the angle is
        measured.

    Returns
    -------
    a: float
        Angle from electrode 1 to electrode 2 in rad, where an angle of 0 is
        parallel to (1,1)->(10,1), an angle of pi measures the angle
        (10,1)->(1,1). The angle is in the range [0..2pi).
    """
    el1_pos = get_elec_positions(electrode1)
    el2_pos = get_elec_positions(electrode2)

    dist_row = (el2_pos[0] - el1_pos[0]) * 400.
    dist_col = (el2_pos[1] - el1_pos[1]) * 400.

    a = np.angle(np.complex(dist_row, dist_col))
    if a < 0.:
        a = 2.0 * np.pi + a

    return a


def get_elec_center_angle(electrode):
    """
    Returns the angle from the center of the array to a given electrode.

    Parameters
    ----------
    electrode: int
        Linear electrode IDs (1-100).

    Returns
    -------
    a: float
        Angle from the center to the electrode in rad, where an angle of pi/4
        rad (45 degrees) is defined from center to electrode 100. The angle is
        in the range [0..2pi).
    """
    el_pos = get_elec_positions(electrode)

    dist_row = (el_pos[0] - 5.5) * 400.
    dist_col = (el_pos[1] - 5.5) * 400.

    a = np.angle(np.complex(dist_row, dist_col))
    if a < 0.:
        a = 2.0 * np.pi + a

    return a


def butterworth(data, low, high, order=3, fs=1000):
    """
    Filters an array using the filtfilt function with a Butterworth filter of
    given order.
    """
    nyq = 0.5 * fs
    lowcut, highcut = low / nyq, high / nyq
    l, h = pysig.butter(order, [lowcut, highcut], btype='band')
    buttered_data = pysig.filtfilt(l, h, data)

    return buttered_data


def applyfilter(ansig, lowcut, highcut, order):
    '''
    Apply a butterworth passband filter to an AnalogSignal.

    Parameters
    ----------
    signal : neo.AnalogSignal
        Signal to filter
    lowcut : float
        Lower cutoff frequency in Hz.
    highcut : float
        Upper cut-off frequency in Hz.
    order : int
        Filter order.

    Returns
    -------
    neo.AnalogSignal
        Contains the filtered signal.
    '''

    # TODO: Integrate into jelephant, use existing filtering method
    return ansig.duplicate_with_new_array(
        butterworth(
            data=ansig.magnitude, low=lowcut, high=highcut,
            order=order, fs=ansig.sampling_rate.rescale(pq.Hz)))


def applyzscore(ansig):
    '''
    Apply a z-score operation to an AnalogSignal.

    This operation subtracts the mean of the signal, and divides by its
    standard deviation.

    Parameters
    ----------
    signal : neo.AnalogSignal
        Signal to z-score.

    Returns
    -------
    neo.AnalogSignal
        Contains the z-scored signal.
    '''

    return ansig.duplicate_with_new_array(
        (ansig.magnitude - np.mean(ansig.magnitude)) / np.std(ansig.magnitude))


def applyhilbert(ansig):
    '''
    Apply a Hilbert transform to an AnalogSignal in order to obtains its
    (complex) analytic signal.

    Parameters
    -----------
    signal : neo.AnalogSignal
        Signal to transform.

    Returns
    -------
    neo.AnalogSignal
        Contains the analytic signal.
    '''

    # To speed up calculation of the Hilbert transform, make sure we change the
    # signal to be of a length that is a power of two. Failure to do so results
    # in computations of certain signal lengths to not finish (or finish in
    # absurd time).
    n_org = len(ansig.magnitude)
    n_opt = int(math.pow(2, math.ceil(math.log(n_org) / math.log(2))))

    # Right-pad signal to desired length using the signal itself
    s = np.hstack((ansig.magnitude, ansig.magnitude[:n_opt - n_org]))

    return ansig.duplicate_with_new_array(
        pysig.hilbert(s, N=n_opt)[:n_org])



def spike_triggered_phase(spiketrains, hilbert_transform, interpolate, extrasignal):
    """
    Calculate the set of spike-triggered phases of an AnalogSignal.

    Parameters
    ----------
    spiketrains : Spiketrain or list of Spiketrain
        Spiketrains on which to trigger hilbert_transform extraction
    hilbert_transform : AnalogSignal or list of AnalogSignal
        AnalogSignal of the complex analytic signal (e.g., returned by the
        elephant.signal_processing.hilbert()). All spike trains are compared to
        this signal, if only one signal is given. Otherwise, length of
        hilbert_transform must match the length of spiketrains.
    interpolate : bool
        If True, the phases and amplitudes of hilbert_transform for spikes
        falling between two samples of signal is interpolated. Otherwise, the
        closest sample of hilbert_transform is used.

    Returns
    -------
    phases : list of arrays
        Spike-triggered phases. Entries in the list correspond to the
        SpikeTrains in spiketrains. Each entry contains an array with the
        spike-triggered angles (in rad) of the signal.
    amp : list of arrays
        Corresponding spike-triggered amplitudes.
    times : list of arrays
        A list of times corresponding to the signal
        Corresponding times (corresponds to the spike times).

    Example
    -------
    Create a 20 Hz oscillatory signal sampled at 1 kHz and a random Poisson
    spike train:

    >>> f_osc = 20. * pq.Hz
    >>> f_sampling = 1 * pq.ms
    >>> tlen = 100 * pq.s
    >>> time_axis = np.arange(
            0, tlen.rescale(pq.s).magnitude,
            f_sampling.rescale(pq.s).magnitude) * pq.s
    >>> analogsignal = AnalogSignal(
            np.sin(2 * np.pi * (f_osc * time_axis).simplified.magnitude),
            units=pq.mV, t_start=0 * pq.ms, sampling_period=f_sampling)
    >>> spiketrain = elephant.spike_train_generation.
            homogeneous_poisson_process(
                50 * pq.Hz, t_start=0.0 * ms, t_stop=tlen.rescale(pq.ms))

    Calculate spike-triggered phases and amplitudes of the oscillation:
    >>> phases, amps, times = elephant.phase_analysis.spike_triggered_phase(
            spiketrain,
            elephant.signal_processing.hilbert(analogsignal),
            interpolate=True)
    """

    # Convert inputs to lists
    if type(spiketrains) is not list:
        spiketrains = [spiketrains]

    if type(hilbert_transform) is not list:
        hilbert_transform = [hilbert_transform]

    if type(extrasignal) is not list:
        extrasignal = [extrasignal]

    # Number of signals
    num_spiketrains = len(spiketrains)
    num_phase = len(hilbert_transform)

    # For each trial, select the first input
    start = [elem.t_start for elem in hilbert_transform]
    stop = [elem.t_stop for elem in hilbert_transform]

    result_phases = [[] for _ in range(num_spiketrains)]
    result_amps = [[] for _ in range(num_spiketrains)]
    result_times = [[] for _ in range(num_spiketrains)]
    result_extra = [[] for _ in range(num_spiketrains)]

    # Step through each signal
    for spiketrain_i, spiketrain in enumerate(spiketrains):
        # Check which hilbert_transform AnalogSignal to look at - if there is
        # only one then all spike trains relate to this one, otherwise the two
        # lists of spike trains and phases are matched up
        if num_phase > 1:
            phase_i = spiketrain_i
        else:
            phase_i = 0

        # Take only spikes which lie directly within the signal segment -
        # ignore spikes sitting on the last sample
        sttimeind = np.where(np.logical_and(
            spiketrain >= start[phase_i], spiketrain < stop[phase_i]))[0]

        # Find index into signal for each spike
        ind_at_spike = np.round(
            (spiketrain[sttimeind] - hilbert_transform[phase_i].t_start) /
            hilbert_transform[phase_i].sampling_period).magnitude.astype(int)

        # Extract times for speed reasons
        times = hilbert_transform[phase_i].times

        # Step through all spikes
        for spike_i, ind_at_spike_j in enumerate(ind_at_spike):
            # Difference vector between actual spike time and sample point,
            # positive if spike time is later than sample point
            dv = spiketrain[sttimeind[spike_i]] - times[ind_at_spike_j]

            # Make sure ind_at_spike is to the left of the spike time
            if dv < 0 and ind_at_spike_j > 0:
                ind_at_spike_j = ind_at_spike_j - 1

            if interpolate:
                # Get relative spike occurrence between the two closest signal
                # sample points
                # if z->0 spike is more to the left sample
                # if z->1 more to the right sample
                z = (spiketrain[sttimeind[spike_i]] - times[ind_at_spike_j]) /\
                    hilbert_transform[phase_i].sampling_period

                # Save hilbert_transform (interpolate on circle)
                p1 = np.angle(hilbert_transform[phase_i][ind_at_spike_j])
                p2 = np.angle(hilbert_transform[phase_i][ind_at_spike_j + 1])
                result_phases[spiketrain_i].append(
                    np.angle(
                        (1 - z) * np.exp(np.complex(0, p1)) +
                        z * np.exp(np.complex(0, p2))))

                # Save amplitude
                result_amps[spiketrain_i].append(
                    (1 - z) * np.abs(
                        hilbert_transform[phase_i][ind_at_spike_j]) +
                    z * np.abs(hilbert_transform[phase_i][ind_at_spike_j + 1]))

                result_extra[spiketrain_i].append(
                    extrasignal[phase_i][ind_at_spike_j])

            else:
                p1 = np.angle(hilbert_transform[phase_i][ind_at_spike_j])
                result_phases[spiketrain_i].append(p1)

                # Save amplitude
                result_amps[spiketrain_i].append(
                    np.abs(hilbert_transform[phase_i][ind_at_spike_j]))

                result_extra[spiketrain_i].append(
                    extrasignal[phase_i][ind_at_spike_j])

            # Save time
            result_times[spiketrain_i].append(spiketrain[sttimeind[spike_i]])

    # Convert outputs to arrays
    for i, entry in enumerate(result_phases):
        result_phases[i] = np.array(entry).flatten()
    for i, entry in enumerate(result_amps):
        result_amps[i] = pq.Quantity(entry, units=entry[0].units).flatten()
    for i, entry in enumerate(result_times):
        result_times[i] = pq.Quantity(entry, units=entry[0].units).flatten()
    for i, entry in enumerate(result_extra):
        result_extra[i] = pq.Quantity(entry, units=entry[0].units).flatten()

    return (result_phases, result_amps, result_times, result_extra)


def calc_good_el(block):
    '''
    Returns a list of all good electrodes. Goodnes is defined as one where the
    annotation 'rejMid' is False (i.e., LFP is ok).

    Parameters
    ----------
    block : neo.Block
        Block containing the data.

    Returns
    -------
    numpy.array
        An array of type boolean, where element i is True if the electrode
        i+1 (linear connector aligned ID) should be used in further
        calculations (i.e., corrected and not rejected).
    '''

    a = np.tile([False], (100))

    # Go through all electrodes
    for i in xrange(100):
        # Get the signal
        asig_list = pick.get_sig_analogsignal(
            block, {'ca_id': [i + 1], 'rejIFC': [False]})

        if len(asig_list) != 0:
            a[i] = True

    return a


def _analyze_block(block):
    '''
    Analyze a Block structure to detect if all signals have the same length,
    sampling and t_start. If all AnalogSignals are similar, the corresponding
    values are returned, otherwise an Exception is thrown.

    Parameters
    -----------
    signal : neo.AnalogSignal
        Signal to transform.

    Returns
    -------
    common_len: int
        Length of each AnalogSignal in number of samples
    common_sf: float
        Sampling frequency of each AnalogSignal
    common_ts: float
        Start time of each AnalogSignal
    '''
    common_len = None
    common_sf = None
    common_ts = None

    asig_list = pick.get_sig_analogsignal(block, {'rejMid': [False]})
    for asig in asig_list:
        len_asig = len(asig)
        sr_asig = asig.sampling_rate
        ts_asig = asig.t_start

        if common_len is None:
            common_len = len_asig
            common_sf = sr_asig
            common_ts = ts_asig
        else:
            if common_len != len_asig:
                raise Exception
            if common_sf != sr_asig:
                raise Exception
            if common_ts != ts_asig:
                raise Exception

    return common_len, common_sf, common_ts


def calc_phases(block):
    '''
    Calculate the phases from a neo block containing the Hilbert transforms as
    analog signals.


    Parameters
    ----------
    block : neo.Block
        Block containing the data, i.e., a block of AnalogSignals that contain
        the complex-valued Hilbert transforms.

    Returns
    -------
    neo.AnalogSignalArray
        Contains the phase signal for each electrode. Signals are ordered
        according to their connector aligned linear ID. Signals which are not
        connected (compare calc_good_el) are zero.
    '''

    common_len, common_sf, common_ts = _analyze_block(block)

    # Create an array
    a = np.zeros((common_len, 100))

    # Go through all electrodes
    for i in xrange(100):
        # Get the signal
        asig_list = pick.get_sig_analogsignal(
            block, {'ca_id': [i + 1], 'rejMid': [False]})

        # If signal was obtained (i.e., not rejected or not connected)
        if len(asig_list) == 1:
            a[:, i] = np.angle(asig_list[0].magnitude)

    return neo.AnalogSignalArray(
        a * pq.dimensionless, sampling_rate=common_sf, t_start=common_ts)


def calc_amps(block):
    '''
    Calculate the amplitudes from a neo block containing the Hilbert transforms
    as analog signals.


    Parameters:
    -----------
    block : neo.Block
        Block containing the data.

    Returns:
    --------
    neo.AnalogSignalArray
        Contains the amplitude signal for each electrode. Signals are ordered
        according to their connector aligned linear ID. Signals which are not
        connected (compare calc_good_el) are zero.
    '''

    common_len, common_sf, common_ts = _analyze_block(block)

    # Create an array
    a = np.zeros((common_len, 100))

    # Go through all electrodes
    for i in xrange(100):
        # Get the signal
        asig_list = pick.get_sig_analogsignal(
            block, {'ca_id': [i + 1], 'rejMid': [False]})

        # If signal was obtained (i.e., not rejected or not connected)
        if len(asig_list) == 1:
            a[:, i] = np.abs(asig_list[0].magnitude)

    return neo.AnalogSignalArray(
        a * pq.dimensionless, sampling_rate=common_sf, t_start=common_ts)


def calc_phgr(sigarray, nextneighbor_distance, good_el):
    '''
    Calculate the phase gradients from an AnalogSignalArray containing the
    phases of 100 electrodes.

    Parameters
    ----------
    sigarray : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase time series (i.e., the output of calc_phases).
    nextneighbor_distance : int
        Number of nearest neighbors to consider in calculating the gradient.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignalArray
        Contains the phase gradient signal for each electrode. Signals are
        ordered according to their connector aligned linear ID. Signals which
        are not connected (compare calc_good_el) are zero. The unit is in
        rad/mm, the direction of the arrow is the average direction of the
        gradient.
    '''

    common_len = sigarray.shape[0]

    # Create an array
    result = np.zeros((common_len, 100), dtype=np.complex)

    # Go through all electrodes
    for i in xrange(100):
        if good_el[i]:
            pdiff_row = []
            pdiff_col = []

            # Get time slice
            sig = sigarray[:, i].magnitude

            # Create gradient along row
            nextneighbors = get_nn_row(i + 1, 2)
            for nextneighbor in nextneighbors:
                # Index of next neighbor
                j = nextneighbor - 1

                # Neighbor exists?
                if good_el[j]:
                    # Get the signal
                    nn_sig = sigarray[:, j].magnitude

                    u = norm_angle(nn_sig - sig)
                    dist = get_elec_distance(i + 1, j + 1)

                    # phase gradient between these two electrodes
                    pdiff_row.append(
                        np.sign(j - i) * u / (dist / 1000.))

            # Create gradient along columns
            nextneighbors = get_nn_col(i + 1, 2)
            for nextneighbor in nextneighbors:
                # Index of next neighbor
                j = nextneighbor - 1

                # Neighbor exists?
                if good_el[j]:
                    # Get the signal
                    nn_sig = sigarray[:, j].magnitude

                    u = norm_angle(nn_sig - sig)
                    dist = get_elec_distance(i + 1, j + 1)

                    # phase gradient between these two electrodes
                    pdiff_col.append(
                        np.sign(j - i) * u / (dist / 1000.))
            if len(pdiff_row) == 0:
                pdiff_row = [np.zeros(common_len)]
            if len(pdiff_col) == 0:
                pdiff_col = [np.zeros(common_len)]

            result[:, i] = np.mean(np.array(pdiff_row), axis=0) + \
                np.complex(0, 1) * np.mean(np.array(pdiff_col), axis=0)

    return neo.AnalogSignalArray(
        result / pq.mm,
        sampling_rate=sigarray.sampling_rate, t_start=sigarray.t_start)


def calc_phgr_stencil(sigarray, good_el):
    '''
    Calculate the phase gradients from an AnalogSignalArray containing the
    phases of 100 electrodes.

    NOTE: This function implements the numerical methods for calculating the
    gradient. However: this only works for non-circular values -- it is unclear
    how these formulas would need to be adapted... Maybe by unwrapping the
    phase maps first, i.e., adding 360 degree if necessary?

    Parameters
    ----------
    sigarray : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase time series (i.e., the output of calc_phases).
    nextneighbor_distance : int
        Number of nearest neighbors to consider in calculating the gradient.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignalArray
        Contains the phase gradient signal for each electrode. Signals are
        ordered according to their connector aligned linear ID. Signals which
        are not connected (compare calc_good_el) are zero. The unit is in
        rad/mm, the direction of the arrow is the average direction of the
        gradient.
    '''

    common_len = sigarray.shape[0]

    # Create an array
    result = np.zeros((common_len, 100), dtype=np.complex)

    # Go through all electrodes
    for i in xrange(100):
        if good_el[i]:
            # Create gradient along row
            nextneighbors = get_nn_row(i + 1, 2)
            if len(nextneighbors) == 5 and \
                    all([good_el[_ - 1] for _ in nextneighbors]):
                # Five-point stencil method
                drow = (
                    sigarray[:, nextneighbors[0] - 1].magnitude -
                    8. * sigarray[:, nextneighbors[1] - 1].magnitude +
                    8. * sigarray[:, nextneighbors[2] - 1].magnitude -
                    sigarray[:, nextneighbors[3] - 1].magnitude) / 12. / .4
            else:
                nextneighbors = get_nn_row(i + 1, 1)
                if len(nextneighbors) == 2 and \
                        all([good_el[_ - 1] for _ in nextneighbors]):
                    # First-order method near edge or bad electrode
                    drow = (
                        sigarray[:, nextneighbors[1] - 1].magnitude -
                        sigarray[:, nextneighbors[0] - 1].magnitude) / 2. / .4
                elif good_el[nextneighbors[0] - 1]:
                    # Poor-man's gradient at edge
                    drow = np.sign(nextneighbors[0] - (i + 1)) * (
                        sigarray[:, nextneighbors[0] - 1].magnitude -
                        sigarray[:, i].magnitude) / .4
                else:
                    drow = 0.

            # Create gradient along column
            nextneighbors = get_nn_col(i + 1, 2)
            if len(nextneighbors) == 5 and \
                    all([good_el[_ - 1] for _ in nextneighbors]):
                # Five-point stencil method
                dcol = (
                    sigarray[:, nextneighbors[0] - 1].magnitude -
                    8. * sigarray[:, nextneighbors[1] - 1].magnitude +
                    8. * sigarray[:, nextneighbors[2] - 1].magnitude -
                    sigarray[:, nextneighbors[3] - 1].magnitude) / 12. / .4
            else:
                nextneighbors = get_nn_col(i + 1, 1)
                if len(nextneighbors) == 2 and \
                        all([good_el[_ - 1] for _ in nextneighbors]):
                    # First-order method near edge or bad electrode
                    dcol = (
                        sigarray[:, nextneighbors[1] - 1].magnitude -
                        sigarray[:, nextneighbors[0] - 1].magnitude) / 2. / .4
                elif good_el[nextneighbors[0] - 1]:
                    # Poor-man's gradient at edge
                    dcol = np.sign(nextneighbors[0] - (i + 1)) * (
                        sigarray[:, nextneighbors[0] - 1].magnitude -
                        sigarray[:, i].magnitude) / .4
                else:
                    dcol = 0.

            result[:, i] = np.complex(norm_angle(drow), norm_angle(dcol))

    return neo.AnalogSignalArray(
        result / pq.mm,
        sampling_rate=sigarray.sampling_rate, t_start=sigarray.t_start)


def calc_phlat(sigarray, reference_phase, good_el):
    '''
    Calculate the phase latency maps (Muller et al, 2014) from an
    AnalogSignalArray containing the phases of 100 electrodes.

    The phase latency method calculates, channel by channel, the time it takes
    for the phase signal on that channel to cross a certain reference phase.
    The phase crossing is evaluated by linearly interpolating time between the
    phase value before and after the crossing. The final resultant map is
    shifted in time such that the lowest latency is zero.

    Parameters
    ----------
    sigarray : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase time series.
    reference_phase : float
        Reference phase crossing up to which to calculate the time difference.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignalArray
        Contains the phase latency maps for each electrode. Signals are
        ordered according to their connector aligned linear ID, excluding those
        signals which are not connected (compare calc_good_el). The unit is in
        ms.
    '''

    common_len = sigarray.shape[0]

    # create an array
    result = np.zeros((common_len, 100))

    # go through all electrodes
    for i in xrange(100):
        if good_el[i]:
            # get time slice
            sig = sigarray[:, i].magnitude

            # find all crossings of the reference phase and get the indices of
            # the bin before each crossing
            cross_idx = np.where(
                np.array(np.diff(sig > reference_phase) == 1, dtype=int))[0]

            # location of last crossing, to indicate which parts of the signal
            # need to be calculated
            last_cross = 0

            # step through all crossings, and calculate the latency up to each
            # crossing for all bins immediately before
            for cross_i in cross_idx:
                # calculate exact offset of crossing by linear interpolation
                perc_cross = (
                    reference_phase - sig[cross_i]) / (
                    sig[cross_i + 1] - sig[cross_i])

                # calculate latency: number of bins before crossing, plus
                # interpolated crossing
                result[last_cross:cross_i + 1, i] = \
                    cross_i + perc_cross - np.arange(last_cross, cross_i + 1)

                last_cross = cross_i + 1

    return neo.AnalogSignalArray(
        result * sigarray.sampling_period.units,
        sampling_rate=sigarray.sampling_rate, t_start=sigarray.t_start)


def calc_cohphgr(sigarray, nextneighbor_distance, good_el):
    '''
    Calculate the coherence phase gradients from an AnalogSignalArray
    containing the phases of 100 electrodes.

    Parameters:
    -----------
    sigarray : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    nextneighbor_distance : int
        Number of nearest neighbors to consider in calculating the gradient.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns:
    --------
    neo.AnalogSignalArray
        Contains the coherence phase gradient signal for each electrode.
        Signals are ordered according to their connector aligned linear ID,
        excluding those signals which are not connected (compare calc_good_el).
        The unit is in rad/m.
    '''
    common_len = sigarray.shape[0]

    # create an array
    result = np.zeros((common_len, 100), dtype=np.complex)

    # go through all electrodes
    for i in xrange(100):
        if good_el[i]:

            nextneighbors = get_nn(i + 1, nextneighbor_distance)
            phgrads = []

            # get the signal on own electrode
            sig = sigarray[:, i].magnitude
            phgrads.append(sig / np.abs(sig))

            for nextneighbor in nextneighbors:
                # index of next neighbor
                j = nextneighbor - 1

                # neighbor exists?
                if good_el[j]:
                    # get the signal
                    sig = sigarray[:, j].magnitude
                    phgrads.append(sig / np.abs(sig))

            result[:, i] = np.mean(phgrads, axis=0)

    return neo.AnalogSignalArray(
        result * pq.dimensionless,
        sampling_rate=sigarray.sampling_rate, t_start=sigarray.t_start)


def calc_meas_var_phase(signal, good_el):
    '''
    Calculate the 'variance phase' measure of each wave frame.

    This measure calculates the variance of the phase across all good
    electrodes.

    Parameters
    ----------
    frames : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the phases.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignal
        Contains the variance of the phase across electrodes.
    '''
    return neo.AnalogSignal(
        1.0 - np.abs(np.mean(np.exp(
            np.complex(0, 1) * signal[:, good_el].magnitude), axis=1)) *
        pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_var_vecfield(signal, good_el):
    '''
    Calculate the 'variance phase gradient' measure of each wave frame.

    This measure calculates the average phase gradient across all good
    electrodes, and returns the length of the average vector as measure of the
    variance of gradients.

    Parameters
    ----------
    frames : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignal
        Contains the variance of the phase gradient across electrodes.
    '''
    return neo.AnalogSignal(
        1.0 - np.abs(np.mean(signal[:, good_el].magnitude, axis=1)) *
        pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_var_norm_vecfield(signal, good_el):
    '''
    Calculate the 'variance of normalized phase gradients' measure of each wave
    frame.

    This measure calculates the average of all normalized phase gradients (to
    length 1, but keeping their direction) across all good electrodes, and
    returns the length of the average vector as measure of the variance of
    gradients.

    Parameters
    ----------
    frames : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignal
        Contains the variance of the normalized phase gradient across
        electrodes.
    '''
    return neo.AnalogSignal(
        1.0 - np.abs(np.mean(
            signal[:, good_el].magnitude /
            np.abs(signal[:, good_el].magnitude), axis=1)) * pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_avg_len_vecfield(signal, good_el):
    '''
    Calculate the 'average length of the coherence phase gradient' measure of
    each wave frame.

    This measure calculates the average length of the coherence phase gradients
    across all good electrodes.

    Parameters
    ----------
    frames : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        coherence phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignal
        Contains the average length of the coherence phase gradient across
        electrodes.
    '''
    return neo.AnalogSignal(
        np.mean(np.abs(signal[:, good_el].magnitude), axis=1) *
        pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_contingency(signal, good_el):
    '''
    Calculate the 'contingency' measure of each wave frame.

    This measure checks for each electrode, to which electrode the phase
    gradient points to. Then it checks if the phase gradient of that electrode
    matches that of the originating electrode, indicating a "contingency"
    from one electrode to the next. This difference is measures as the cosine
    distance between phase gradients. The measure returns the average of this
    distance measure across all electrodes.

    Parameters
    ----------
    frames : neo.core.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns:
    --------
    neo.core.AnalogSignal
        Contains the contingeny measure (-1...1).
    '''
    result = np.zeros((signal.shape[0]))

    # electrode ID offsets for different directions, starting at phi=0 (along x
    # axis)
    offsets = [1, 11, 10, 9, -1, -11, -10, -9]

    # consider only electrodes which are not on the border
    cons_el = []
    for i in xrange(1, 9):
        for j in xrange(1, 9):
            cons_el.append(i + j * 10)

    # pre-calculate angles
    angs = np.angle(signal.magnitude)

    # go through all time points
    for t in range(signal.shape[0]):
        contingencies = []

        # step through all electrodes
        for i in cons_el:
            # if this electrode is good
            if good_el[i]:
                # find at time i, where the arrow points to
                # divide the circle in 8 parts shifted by 2pi/16=pi/8
                # then pick n-th eighth of the circle (=8/2pi=4/pi)
                k = np.mod(angs[t, i] + np.pi / 8.0, 2.0 * np.pi)
                n = int(k * 4.0 / np.pi)

                # find square the gradient points to
                offset = offsets[n]

                # is electrode pointed to good?
                if good_el[i + offset]:
                    # calculate how contingent the phase gradients are
                    contingencies.append(np.cos(
                        angs[t, i] - angs[t, i + offset]))
        result[t] = np.mean(contingencies)

    return neo.AnalogSignal(
        result * pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_outwardishness(signal, good_el):
    '''
    Calculate the 'outwardishness' measure of each wave frame.

    This measure calculates for each electrode of the array, whether its phase
    gradient points outward from the array center. The measure for this is the
    distance between the phase gradient and a vector from the array center to
    the corresponding electrode. The outwardishness measure is the average of
    this distance across the array.

    Parameters
    ----------
    frames : neo.core.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.core.AnalogSignal
        Contains the outwardishness  (-1...1).
    '''
    result = np.zeros((signal.shape[0]))

    # pre-calculate angles
    angs = np.angle(signal.magnitude)

    # go through all time points
    for t in range(signal.shape[0]):
        outwardishnesses = []

        # step through all electrodes
        for i in range(signal.shape[1]):
            # if this electrode is good
            if good_el[i]:
                outwardishnesses.append(np.cos(
                    angs[t, i] - get_elec_center_angle(i)))

        result[t] = np.mean(outwardishnesses)

    return neo.AnalogSignal(
        result * pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_circlishness(signal, good_el):
    '''
    Calculate the 'circlishness' measure of each wave frame.

    This measure calculates for each electrode of the array, whether its phase
    gradient points outward from the array center. The measure for this is the
    distance between the phase gradient and a vector from the array center to
    the corresponding electrode. The circlishness measure is the average of
    this distance across the array.

    Parameters
    ----------
    frames : neo.core.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.core.AnalogSignal
        Contains the circlishness  (-1...1).
    '''
    result = np.zeros((signal.shape[0]))

    # pre-calculate angles
    angs = np.angle(signal.magnitude)

    # go through all time points
    for t in range(signal.shape[0]):
        circlishnesses = []

        # step through all electrodes
        for i in range(signal.shape[1]):
            # if this electrode is good
            if good_el[i]:
                circlishnesses.append(np.sin(
                    angs[t, i] - get_elec_center_angle(i)))

        result[t] = np.mean(circlishnesses)

    return neo.AnalogSignal(
        result * pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_avg_amp(signal, good_el):
    '''
    Calculate the 'average amplitude' measure of each wave frame.

    This measure calculates the average amplitude (i.e., the envelope of the
    analytic signal) across all good electrodes.

    Parameters
    ----------
    frames : neo.AnalogSignalArray
        Signal array containing the amplitudes for each of the 100 electrodes.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignal
        Contains the average amplitude across electrodes.
    '''
    return neo.AnalogSignal(
        np.mean(signal[:, good_el], axis=1) * signal.units,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_avg_velocity(signal, good_el, f):
    '''
    Calculate the 'average velocity' measure of each wave frame.

    This measure calculates the average wave velocity from the phase gradients
    of all good electrodes.

    Parameters
    ----------
    frames : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).
    f : quantities.Quantity
        Frequency of oscillation of the wave.

    Returns
    -------
    neo.AnalogSignal
        Contains the average velocity of the wave across electrodes (in cm/s).
    '''
    return neo.AnalogSignal(
        2. * np.pi * f /
        np.mean(np.abs(signal[:, good_el].rescale(1. / pq.cm).magnitude),
                axis=1) * pq.cm,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def calc_meas_avg_velocity_directed(signal_phase, signal_phgr, good_el, f):
    '''
    Calculate the 'average velocity directed' measure of each wave frame.

    This measure calculates the wave velocity from the phases and phase
    gradients of all good electrodes, taking into account only the velocity in
    direction of the phase gradient.

    Parameters
    ----------
    signal_phase : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phases.
    signal_phgr : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).
    f : quantities.Quantity
        Frequency of oscillation of the wave.

    Returns
    -------
    neo.AnalogSignal
        Contains the velocity of the wave in direction of the phase gradient
        across electrodes (in cm/s).
    '''
    result = np.zeros((signal_phase.shape[0]))

    # electrode ID offsets for different directions, starting at phi=0 (along x
    # axis)
    offsets = [1, 11, 10, 9, -1, -11, -10, -9]

    # consider only electrodes which are not on the border
    cons_el = []
    for i in xrange(1, 9):
        for j in xrange(1, 9):
            cons_el.append(i + j * 10)

    # pre-calculate angles
    angs = np.angle(signal_phgr.magnitude)
    phases = signal_phase.magnitude

    # go through all time points
    for t in range(signal_phase.shape[0]):
        velocities = []

        # step through all electrodes
        for i in cons_el:
            # if this electrode is good
            if good_el[i]:
                # find at time i, where the arrow points to
                # divide the circle in 8 parts shifted by 2pi/16=pi/8
                # then pick n-th eighth of the circle (=8/2pi=4/pi)
                k = np.mod(angs[t, i] + np.pi / 8.0, 2.0 * np.pi)
                n = int(k * 4.0 / np.pi)

                # find square the gradient points to
                j = i + offsets[n]

                # is electrode pointed to good?
                if good_el[j]:
                    # u=phase difference
                    # dist = distance in um
                    u = norm_angle(phases[t, j] - phases[t, i])
                    dist = get_elec_distance(i + 1, j + 1)

                    # calculate how contingent the phase gradients are
                    # in rad/cm
                    velocities.append(u / (dist / 10000.))
        result[t] = np.mean(velocities)

    return neo.AnalogSignal(
        2. * np.pi * f * result * pq.cm,
        tstart=signal_phase.t_start, sampling_rate=signal_phase.sampling_rate)


def calc_meas_avg_dir(signal, good_el):
    '''
    Calculate the 'average wave direction' measure of each wave frame.

    This measure calculates the average direction of a (planar) wave from the
    phase gradient across all good electrodes. Here, the direction is
    calculated from normalized gradients, i.e., each electrode contributes its
    direction with equal weight, independent of the magnitude of the gradient.

    Parameters
    ----------
    frames : neo.AnalogSignalArray
        Signal array containing for each of the 100 electrodes the (complex)
        phase gradients.
    good_el : numpy.array
        Array of 100 boolean values indicating which of the 100 electrodes is
        good (in connector aligned order).

    Returns
    -------
    neo.AnalogSignal
        Contains the average direction of the wave (in rad).
    '''
    return neo.AnalogSignal(
        np.angle(np.mean(
            signal[:, good_el].magnitude /
            np.abs(signal[:, good_el].magnitude), axis=1)) * pq.dimensionless,
        tstart=signal.t_start, sampling_rate=signal.sampling_rate)


def classify_wave_patterns_thresh(meas, thresh):
    '''
    Classify the pattern of wave propagation using thresholds.

    Parameters
    ----------
    meas : dict of neo.AnalogSignal
        The measures to use in the clustering.

    Returns
    -------
    neo.AnalogSignal
        Sequence of clusters (0-5):
            0: unclassified
            1: planar
            2: synchronized
            3: unstructured
            4: circular
            5: radial
    '''
    # conditions 1-8
    c1 = meas['meas_var_norm_phgr'].magnitude < thresh[1][0]
    c2 = meas['meas_var_norm_phgr'].magnitude >= thresh[1][1]
    c3 = meas['meas_var_phase'].magnitude < thresh[0][0]
    c4 = meas['meas_var_phase'].magnitude >= thresh[0][1]
    c5 = meas['meas_avg_cohphgr'].magnitude < thresh[2][0]
    c6 = meas['meas_contingency_cohphgr'].magnitude >= thresh[3][0]
    c7 = np.abs(meas['meas_circlishness_cohphgr'].magnitude) >= thresh[4][0]
    c8 = np.abs(meas['meas_outwardishness_cohphgr'].magnitude) >= thresh[5][0]

    # state vector (determine length from any of the measures)
    T = np.zeros(len(meas['meas_var_norm_phgr'].magnitude))

    # planar wave if phase gradients are clearly in one direction
    T[c1] = 1

    # otherwise, outward if high degree of outwardishness (inside or outside)
    T[np.logical_and(T == 0, c8)] = 5

    # if not planar wave, detect synchronized state if phases are all the same
    # (c3), and not planar(c2)
    T[np.logical_and(T == 0, np.logical_and(c3, c2))] = 2

    # spiral, if random phases with high contingency and orthogonal
    # outwardishness
    T[np.logical_and(T == 0, np.logical_and(np.logical_and(np.logical_and(
        c2, c4), c6), c7))] = 4

    # random, if phases are all different (c4), not planar (c2), locally
    # not clear in direction (c5)
    T[np.logical_and(T == 0, np.logical_and(np.logical_and(
        c2, c4), c5))] = 3

    return neo.AnalogSignal(T * pq.dimensionless,
                            tstart=meas['meas_var_norm_phgr'].t_start,
                            sampling_rate=meas['meas_var_norm_phgr'].sampling_rate)


def classify_wave_patterns_km(meas):
    '''
    Classify the pattern of wave propagation using kmeans

    Parameters
    ----------
    meas : dict of neo.AnalogSignal
        The measures to use in the clustering.

    Returns
    -------
    neo.AnalogSignal
        Clustered sequence of clusters (1-5).
    neo.AnalogSignal
        Array containing the coordinates of the centroids.
    '''

    # construct an input vector of all measures
    X = np.vstack([meas[i].magnitude for i in [
        'meas_var_phase',
        'meas_var_norm_phgr',
        'meas_avg_cohphgr',
        # TODO!!! Should be phgr!
        #'meas_contingency_phgr',
        'meas_contingency_cohphgr',
        'meas_outwardishness_cohphgr']]).T

    # number of clusters
    cl_no = 5

    whitened_X = vqcluster.whiten(X)
    C, T = vqcluster.kmeans2(whitened_X, cl_no)

    return (
        neo.AnalogSignal(
            T * pq.dimensionless,
            tstart=meas['meas_var_phase'].t_start,
            sampling_rate=meas['meas_var_phase'].sampling_rate),
        C)


def classify_wave_patterns_hc(meas, ds):
    '''
    Classify the pattern of wave propagation using hierarchical clustering.

    Parameters
    ----------
    meas : dict of neo.AnalogSignal
        The measures to use in the clustering.
    ds : int
        Integer indicating the downsampling factor.

    Returns
    -------
    neo.AnalogSignal
        Clustered sequence of clusters (1-5).
    neo.AnalogSignal
        Linkage (ward).
    float :
        Distance between the last two clusters (i.e., all hierarchical clusters
        with larger distances are ignored)
    '''

    # construct an input vector of all measures
    X = np.vstack([meas[i].magnitude[0:-1:ds] for i in [
        'meas_var_phase',
        'meas_var_norm_phgr',
        'meas_avg_cohphgr',
        # TODO
        #'meas_contingency_phgr',
        'meas_contingency_cohphgr',
        'meas_outwardishness_cohphgr']]).T

    # number of clusters
    cl_no = 5

    # calculate ward distances
    Z = hcluster.ward(X)

    # cluster the distance tree
    T = hcluster.fcluster(Z, cl_no, 'maxclust')

    # distance between cluster 4 and 5?
    ct = Z[-(cl_no - 1), 2]

    return (
        neo.AnalogSignal(
            T * pq.dimensionless,
            tstart=meas['meas_var_phase'].t_start,
            sampling_rate=meas['meas_var_phase'].sampling_rate / ds),
        neo.AnalogSignal(
            Z * pq.dimensionless,
            tstart=meas['meas_var_phase'].t_start,
            sampling_rate=meas['meas_var_phase'].sampling_rate / ds),
        ct)


def create_class_epochs(clusters, mindur=1, min_inter_state=20):
    '''
    Create EpochArray from state vector

    Parameters
    ----------
    clusters : neo.AnalogSignal
        State vector with states 0-5, where state 0 is the unclassified
        (ignored) state.
    mindur : int
        The minimum number of consecutive states are required to add
        an Epoch.
    min_inter_state_dur: int
        The amount of time points that the end of one state and the beginning
        of the next state may be separated such that they are added to the
        transition matrix.

    Returns
    -------
    state_epochs : list of neo.EpochArray
        List of 0 EpochArrays for state 0-5
    neo.AnalogSignal
        New statevector with states 0-5, but keeping only those states that
        consist of at least mindur consecutive entries.
    transitions : np.array
        trans[i,j] indicated the numbers of transitions array from state i to
        state j. Probabilities can be obtained by dividing every row i by the
        length of state_epochs[i].times.
    '''

    # times and durations of the 5 Eopchs
    eps_t = {}
    eps_d = {}

    # updated state vector to output
    clusterout = np.zeros((len(clusters.magnitude)))

    # step through all 5 clusters and create empty epochs
    for c in range(6):
        eps_t[c] = np.array([])
        eps_d[c] = np.array([])

    # transition matrix
    transmat = np.zeros((6, 6))

#     # get differences
#     diffs = np.diff(clusters.magnitude)
#
#     # these are the indexes,where state changes occur.
#     # add 1 to get the first index of the new state, not the last index of the
#     # previous state
#     positions = np.where(diffs != 0)[0] + 1
#
#     # get lengths of epochs in steps between positions
#     epochlen = np.diff(positions)

    # get only those positions, whose following epoch is of minimum duration,
    # and which belong to a cluster > 0
#     finalepochs = positions[
#         np.where(
#             np.logical_and(
#                 epochlen > mindur,
#                 clusters.magnitude[positions[0:-1]] > 0))]
#
#     for i in finalepochs:
#         currentclus = clusters.magnitude[i]
#
#         eps_t[currentclus - 1] = np.append(
#             eps_t[currentclus - 1],
#             clusters.times[i].magnitude)
#
#         eps_d[currentclus - 1] = np.append(
#             eps_d[currentclus - 1],
#             clusters.times[i].magnitude - clusters.times[i].magnitude)

    # saves the last detected state of minimum length in order to make the
    # entry into the transition matrix
    last_state = -1
    last_state_end = 0

    # step through all time points
    t1 = 0

    while t1 < len(clusters.magnitude) - 1:
        currentclus = int(clusters.magnitude[t1])
        t2 = t1 + 1
        while t2 < len(clusters) - 1 and clusters.magnitude[t2] == currentclus:
            t2 += 1

        # is the state long enough?
        if (t2 - t1) >= mindur:
            eps_t[currentclus] = np.append(
                eps_t[currentclus],
                clusters.times[t1].magnitude)

            eps_d[currentclus] = np.append(
                eps_d[currentclus],
                clusters.times[t2].magnitude - clusters.times[t1].magnitude)

            # enter this state into the new state vector
            clusterout[t1:t2] = currentclus

            if currentclus > 0:
                # enter state into the transition matrix
                if last_state != -1 and \
                        (t2 - last_state_end) <= min_inter_state:
                    transmat[last_state, currentclus] += 1

                # this state becomes the last state
                last_state = currentclus
                last_state_end = t2

        t1 = t2

    # normalize transition matrix to probability
#     for i in range(6):
#         for j in range(6):
#             transmat[i, j] = 1. * transmat[i, j] / eps_d[i]

    # create EpochArrays
    epochlist = []
    for c in range(6):
        epochlist.append(
            neo.EpochArray(
                times=eps_t[c] * clusters.times.units,
                durations=eps_d[c] * clusters.times.units))

    return (
        epochlist,
        neo.AnalogSignal(
            clusterout * pq.dimensionless,
            tstart=clusters.t_start, sampling_rate=clusters.sampling_rate),
        transmat)


def cut_asa(ep, fp, reset_times, i):
    """
    Helper function to cut an AnalogSignalArray object

    Parameter
    ---------
    ep : neo.EpochArray
    fp : neo.AnalogSignal
    reset_times : bool
    i : int
        Actual index of the loop

    Returns
    -------
    neo.AnalogSignalArray:
        A cut AnalogSignalArray object.
    None:
        Only if the trial doesn't match inside the range and thus a
        AnalogSignalArray object couldn't be constructed.
    """
    # Range between edges
    ind1 = int((
        (ep.times[i] - fp.t_start) /
        (fp.t_stop - fp.t_start)).rescale('dimensionless') * fp.shape[0])
    ind2 = int((
        ind1 + ep.durations[i].rescale(fp.t_start.units) /
        fp.sampling_period).base)
    # Don't trespass edge
    if ind2 > fp.t_start.magnitude + fp.shape[0]:
        ind2 = fp.t_start.magnitude + fp.shape[0]

    # setting zero time to beginning of analogsignal
    if type(reset_times) == bool and reset_times is True:
        (t_start, t_stop) = (0 * ep.durations[i].units, ep.durations[i])

    # keep time stamps of original spiketrain
    elif type(reset_times == bool) and reset_times is False:
        (t_start, t_stop) = (ep.times[i], ep.times[i] + ep.durations[i])
    else:
        raise TypeError("Provided data type %s for reset times is "
                        "not supported" % (type(reset_times)))
    analog_sig = neo.AnalogSignalArray(fp[ind1:ind2, :],
                                       t_stop=t_stop,
                                       sampling_rate=fp.sampling_rate,
                                       units=fp.units,
                                       t_start=t_start,
                                       trial_id=ep.annotations['trial_id'][i],
                                       **fp.annotations)
    return analog_sig


def cut_as(ep, fp, reset_times, i):
    """
    Helper function to cut an AnalogSignal object

    Parameter
    ---------
    ep : neo.EpochArray
    fp : neo.AnalogSignal
    reset_times : bool
    i : int
        Actual index of the loop

    Returns
    -------
    neo.AnalogSignal:
        A cut AnalogSignal object.
    None:
        Only if the trial doesn't match inside the range and thus a
        AnalogSignal object couldn't be constructed.
    """
    # Range between edges
    ind1 = int((
        (ep.times[i] - fp.t_start) /
        (fp.t_stop - fp.t_start)).rescale('dimensionless') * fp.shape[0])
    ind2 = int((
        ind1 + ep.durations[i].rescale(fp.t_start.units) /
        fp.sampling_period).base)
    # Don't trespass edge
    if ind2 > fp.t_start.magnitude + fp.shape[0]:
        ind2 = fp.t_start.magnitude + fp.shape[0]

    # setting zero time to beginning of analogsignal
    if type(reset_times) == bool and reset_times is True:
        (t_start, t_stop) = (0 * ep.durations[i].units, ep.durations[i])

    # keep time stamps of original spiketrain
    elif type(reset_times == bool) and reset_times is False:
        (t_start, t_stop) = (ep.times[i], ep.times[i] + ep.durations[i])
    else:
        raise TypeError("Provided data type %s for reset times is "
                        "not supported" % (type(reset_times)))
    analog_sig = neo.AnalogSignal(fp[ind1:ind2],
                                  t_stop=t_stop,
                                  sampling_rate=fp.sampling_rate,
                                  units=fp.units,
                                  t_start=t_start,
                                  trial_id=ep.annotations['trial_id'][i],
                                  **fp.annotations)
    return analog_sig
