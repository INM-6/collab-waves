"""
Calculate the average, min, max and std deviation of the time points of trial
events

Output files:
    <session>_<filter>_events_param.h5:
        Contains all analysis parameters used.

    <session>_<filter>_events.h5:
        hdf5 file with calculated results
            events:
                dictionary containing all events of a specific monkey,
                condition, and separated by SG
            event_stat:
                dictionary containing lists with the corresponding mean, std,
                min and max
"""
import os

import numpy as np

import quantities as pq

from rg.rgio import ReachGraspIO

import common.projctrl as projctrl
import common.h5py_wrapper.wrapper as h5pyw

import common.wave_main as wave_main
import common.wave_plots as wave_plots

__updated__ = "2018-02-19"


def calc_events(job_id, selected_subsession, selected_filter):
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

    # =========================================================================
    # Set parameters
    # =========================================================================

    selected_subsession_name = selected_subsession[0]
    selected_subsession_cnd = selected_subsession[1]
    selected_monkey = wave_main.monkey_names[selected_subsession[0][0]]

    print(
        "Started job %i: %s - %s" %
        (job_id, selected_subsession[0], selected_filter))

    # Get directories where files are stored
    data_dir, metadata_dir = wave_main.get_data_dir(selected_subsession_name)

    # =========================================================================
    # Load data
    # =========================================================================

    print("Loading data - %s" % selected_subsession_name)
    subsessionobj = ReachGraspIO(
        os.path.join(data_dir, selected_subsession_name),
        metadata_dir=metadata_dir, print_diagnostic=True)

    # Read in data (lazy, all we need are trials times)
    neo_block = subsessionobj.read_block(
        n_starts=[None], n_stops=[None], channel_list=range(1, 97),
        nsx=None, units=None, events=True)

    # Get good electrodes
    results['good_el_id'] = np.where(wave_main.calc_good_el(neo_block))[0] + 1

    for s in wave_plots.event_names:
        if s != 'TS':
            event_trans = {
                'WS': 'FP-ON', 'CUE-ON': 'CUE-ON', 'CUE-OFF': 'CUE-OFF',
                'GO': 'GO-ON', 'SR': 'SR', 'RW': 'RW'}
            s_translate = event_trans[s]
            subsessionobj.add_trials_between(
                neo_block, 'trigger_epochs',
                'TS', s_translate, aligned_trigger='TS')

            ep = neo_block.segments[0].epocharrays[0]

            trial_types = ep.annotations['trial_type']
            ts = ep.durations.rescale(pq.ms).magnitude

            results['events'][selected_monkey][
                selected_subsession_cnd]['all'][s].extend(list(ts))
            results['events'][selected_monkey][selected_subsession_cnd]['PG'][
                s].extend(
                    list(ts[np.where(np.in1d(
                        trial_types, [85, 86, 89, 101, 149]))]))
            results['events'][selected_monkey][selected_subsession_cnd]['SG'][
                s].extend(
                    list(ts[np.where(np.in1d(
                        trial_types, [106, 154, 166, 169, 170]))]))

            neo_block.segments[0].epocharrays = []


# =============================================================================
# Main function
# =============================================================================

if __name__ == '__main__':
    # Set up project, do not delete existing data
    pc = projctrl.ProjectControl(
        project_name='reachgrasp-spikewave', script_name='calc_events',
        clear_data=True)

    param = {}

    # Create empty results vector
    results = {}
    results['events'] = {'Tanya': {}, 'Nikos2': {}, 'Lilou': {}}
    for k in results['events'].keys():
        results['events'][k] = {1: {}, 2: {}}
        for l in results['events'][k].keys():
            for j in ['all', 'SG', 'PG']:
                results['events'][k][l][j] = {}
                for m in wave_plots.event_names:
                    results['events'][k][l][j][m] = []

    results['event_stat'] = {'Tanya': {}, 'Nikos2': {}, 'Lilou': {}}
    for k in results['event_stat'].keys():
        results['event_stat'][k] = {1: {}, 2: {}}
        for l in results['event_stat'][k].keys():
            for j in ['all', 'SG', 'PG']:
                results['event_stat'][k][l][j] = {}
                for m in wave_plots.event_names:
                    results['event_stat'][k][l][j][m] = []

    # Fetch job ID
    job_id = wave_main.get_jobid()

    # Create job(s) and prefetch files if required
    jobparameters = wave_main.create_file_filter_task(
        wave_main.selected_datasets, wave_main.selected_filter_names[0:1],
        fetch=False)

    # Launch jobs
    wave_main.run_task(
        calc_events, job_id, jobparameters, spare_core=-1)

    # Calculate statistics
    for k in results['event_stat'].keys():
        for l in results['event_stat'][k].keys():
            for j in ['all', 'PG', 'SG']:
                for m in wave_plots.event_names:
                    print(k, l, j, m, results['events'][k][l][j][m])
                    if results['events'][k][l][j][m]:
                        results['event_stat'][k][l][j][m] = [
                            np.mean(results['events'][k][l][j][m]),
                            np.std(results['events'][k][l][j][m]),
                            np.min(results['events'][k][l][j][m]),
                            np.max(results['events'][k][l][j][m])]
                        print results['event_stat'][k][l][j][m]
                    else:
                        results['event_stat'][k][l][j][m] = [
                            0, 0, 0, 0]

    # Write _results to disk
    print("Saving results file")
    h5pyw.add_to_h5(
        pc.result_path + 'events.h5',
        results, write_mode='w', overwrite_dataset=True)

    # Write parameters to disk
    print("Saving parameters")
    h5pyw.add_to_h5(
        pc.result_path + 'events_param.h5',
        param, write_mode='w', overwrite_dataset=True)
