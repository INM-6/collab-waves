"""
Loads a session and calculates/saves basic properties of the wave-like beta
activity and saves the results. These properties are all saved as 10x10 maps of
the activity.

Output files:
    <session>_<filter>_param.h5:
        Contains all analysis parameters used.

    <session>_<filter>_frames.h5:
        Neo hdf5 file with calculated quantities as Neo AnalogSignalArrays:
            phase_frames: phases
            amp_frames: filtered, z-scored amplitudes
            phgr_frames: phases gradients
            cohphgr_frames: coherence phase gradients

        For all AnalogSignalArrays, the signal number X in  [:,X] is the linear
        electrode ID (bottom left is 0, top right is 99).
        Unconnected electrodes are set to 0.

    <session>_<filter>_info.h5:
        Contains paths for individual data objects into _frames file. In
        addition:
            good_el: Boolean array containing all good electrodes in each frame
"""
import os
import sys

# Add Python path and set matplotlib backend
print("\n".join(sys.path))
sys.path.insert(0, "common")
sys.path.insert(0, "reachgrasp-spikewave/src/common")
import matplotlib
matplotlib.use("Agg")

import numpy as np

import neo.io.hdf5io as nh5
from rgio import ReachGraspIO

import pick as pick
import projctrl as projctrl

import h5py_wrapper.wrapper as h5pyw

import wave_main as wave_main

__updated__ = "2018-04-25"


def calc_waveproperties(job_id, selected_subsession, selected_filter):
    """
    Main function to perform the analysis on one subsession.

    Parameters
    ----------
        job_id: integer
            Numeric job ID
        selected_subsession: string
            Filename of the selected subsession.
        selected_filter: string
            Parameter that specifies the filtering to be used.
    """

    # =============================================================================
    # Set and save parameters
    # =============================================================================
    selected_subsession_name = selected_subsession[0]
    selected_subsession_cnd = selected_subsession[1]
    selected_monkey = wave_main.monkey_names[selected_subsession[0][0]]

    print(
        "Started job %i: %s - %s" %
        (job_id, selected_subsession[0], selected_filter))

    # Get directories where files are stored
    data_dir, metadata_dir = wave_main.get_data_dir(selected_subsession_name)

    param = {}

    param['subsession'] = selected_subsession_name
    param['condition'] = selected_subsession_cnd
    param['monkey'] = selected_monkey

    param['filter'] = selected_filter

    # Extract frequencies from selected filter
    if 'lowcut' in wave_main.selected_filters[selected_filter]:
        # Generic filter band
        param['lowcut'] = wave_main.selected_filters[
            selected_filter]['lowcut']
        param['highcut'] = wave_main.selected_filters[
            selected_filter]['highcut']
        param['order'] = wave_main.selected_filters[
            selected_filter]['order']
    else:
        # Monkey specific filter band
        param['lowcut'] = wave_main.selected_filters[selected_filter][
            selected_monkey]['lowcut']
        param['highcut'] = wave_main.selected_filters[selected_filter][
            selected_monkey]['highcut']
        param['order'] = wave_main.selected_filters[selected_filter][
            selected_monkey]['order']
    param['frequency'] = (param['lowcut'] + param['highcut']) / 2.

    param['nextneighbor_distance'] = 2
    param['alignment'] = 'Connector'

    # =============================================================================
    # Load data
    # =============================================================================

    print("Loading data - %s" % selected_subsession_name)
    subsessionobj = ReachGraspIO(
        os.path.join(data_dir, selected_subsession_name),
        metadata_dir=metadata_dir, print_diagnostic=True)

    # Read in data
    neo_block = subsessionobj.read_block(
        n_starts=[None], n_stops=[None], channel_list=range(1, 97),
        nsx=2, events=True)

    # =============================================================================
    # Process data
    # =============================================================================

    # Save various results and path location of frames in neo hdf5
    results = {}

    # Open Neo hdf5 file for saving results
    framesfile = pc.result_path + selected_subsession_name + '_filter_' + \
        selected_filter + '_neo_frames.h5'
    if os.path.exists(framesfile):
        os.remove(framesfile)
    ho = nh5.NeoHdf5IO(filename=framesfile)

    print("Filtering blocks - %s" % selected_subsession_name)
    pick.map_as(
        neo_block, wave_main.applyfilter,
        annotations=None, lowcut=param['lowcut'],
        highcut=param['highcut'], order=param['order'])

    print("z-scoring block - %s" % selected_subsession_name)
    pick.map_as(
        neo_block, wave_main.applyzscore, annotations=None)

    print("Calculating the analytic  - %s" % selected_subsession_name)
    pick.map_as(
        neo_block, wave_main.applyhilbert, annotations=None)

    print("Determining connected, non-broken, non-rejected electrodes - %s" %
          selected_subsession_name)

    # This is a micro solution to obtain the broken electrodes
    if selected_monkey in ['Tanya', 'Nikos2', 'Lilou']:
        import h5py_wrapper.wrapper as h5py
        try:
            rej = h5py.load_h5(os.path.join(
                metadata_dir, "source", "rejections", "data",
                selected_subsession_name + "_lfp_el_rej_015Hz-030Hz-03o.hdf5"))
            list_of_rej_el = rej['bad_electrodes']
        except:
            rej = h5py.load_h5(os.path.join(
                metadata_dir, "source", "rejections", "data",
                selected_subsession_name + "_rej_012Hz-040Hz-03o.hdf5"))
            list_of_rej_el = rej['elrej']['bad_electrodes']
        print(list_of_rej_el)
        for asig in neo_block.segments[0].analogsignals:
            if asig.annotations['electrode_id'] in list_of_rej_el:
                asig.annotations['rejIFC'] = True
            else:
                asig.annotations['rejIFC'] = False
    if selected_monkey == 'Lilou':
        list_of_rej_el = [33, 83, 92]
        for asig in neo_block.segments[0].analogsignals:
            if asig.annotations['ca_id'] in list_of_rej_el:
                asig.annotations['rejIFC'] = True
            else:
                asig.annotations['rejIFC'] = False

    results['good_el'] = wave_main.calc_good_el(neo_block)

    # Check number of non-broken horizontal/vertical neighbors is large,
    # required to calculate meaningful gradients
    for i in range(1, 101):
        nns = wave_main.get_nn_nodiag(i, param['nextneighbor_distance'])
        if np.count_nonzero(
                np.array([results['good_el'][_ - 1] for _ in nns])) < 3 and \
                results['good_el'][i - 1]:
            results['good_el'][i - 1] = False
            # Restart search, since we marked one electrode as bad
            i = 1

    print("Calculating phase  - %s" % selected_subsession_name)
    phase_frames = wave_main.calc_phases(neo_block)
    ho.save(phase_frames)
    results['phase_frames_path'] = phase_frames.hdf5_path

    print("Calculating amplitude frames - %s" % selected_subsession_name)
    amp_frames = wave_main.calc_amps(neo_block)
    ho.save(amp_frames)
    results['amp_frames_path'] = amp_frames.hdf5_path

    print("Calculating phase gradient frames - %s" % selected_subsession_name)
    phgr_frames = wave_main.calc_phgr(
        phase_frames, param['nextneighbor_distance'], results['good_el'])
    ho.save(phgr_frames)
    results['phgr_frames_path'] = phgr_frames.hdf5_path

    print("Calculating coherence phase gradient frames - %s" %
          selected_subsession_name)
    cohphgr_frames = wave_main.calc_cohphgr(
        phgr_frames, param['nextneighbor_distance'], results['good_el'])
    ho.save(cohphgr_frames)
    results['cohphgr_frames_path'] = cohphgr_frames.hdf5_path

    ho.close()

    # Write file containing the table of contents of the Neo hdf5
    print("Saving info file - %s" % selected_subsession_name)
    h5pyw.add_to_h5(
        pc.result_path + selected_subsession_name +
        '_filter_' + selected_filter +
        '_info.h5', results, write_mode='w', overwrite_dataset=True)

    # Write parameters to disk
    print("Saving parameters - %s" % selected_subsession_name)
    h5pyw.add_to_h5(
        pc.result_path + selected_subsession_name +
        '_filter_' + selected_filter +
        '_param.h5', param, write_mode='w', overwrite_dataset=True)


# =============================================================================
# Main function
# =============================================================================

if __name__ == '__main__':
    # Set up project, do not delete existing data
    pc = projctrl.ProjectControl(
        project_name='collab-waves', 
        script_name='calc_waveproperties',
        clear_data=False)

    # Fetch job ID
    job_id = wave_main.get_jobid()

    # Create job(s) and prefetch files if required
    jobparameters = wave_main.create_file_filter_task(
        wave_main.selected_datasets, wave_main.selected_filter_names[0:1],
        fetch=False)

    # Launch jobs (suggest: spare_core=6)
    wave_main.run_task(
        calc_waveproperties, job_id, jobparameters, spare_core=6)
