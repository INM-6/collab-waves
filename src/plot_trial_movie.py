"""
Plots a movie of a given trial.
"""

import os
import colorsys

import numpy as np
import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.image
import matplotlib.offsetbox

import quantities as pq
import neo.io.hdf5io as nh5

# Elephant to perform phase extraction
import elephant.signal_processing

import pick as pick
from rgio import ReachGraspIO

import projctrl as projctrl
import h5py_wrapper.wrapper as h5pyw

import wave_main as wave_main
import wave_plots as wave_plots

__updated__ = "2018-04-27"


# =============================================================================
# Frame function
# =============================================================================

def phasemovie_frames(frame_i, t_i, frame_interval):
    """
    Plot the frame number frame_i, time t_i, at a frame interval of
    frame_interval.
    """

    # Find index into signals at time t_i
    frame_idx = np.argmin(np.abs((t_i - signal_cut.times).magnitude))

    # Offset for plotting, since all signals start at time zero
    tpre_ms = param['tpre'].rescale(pq.ms)
    tpre = tpre_ms.magnitude

    # =========================================================================
    # Title
    # =========================================================================

    # TODO: Remove if this complete block if it always works
    x = np.where(
        neo_block.segments[0].epocharrays[0].annotations['trial_id'] ==
        param['trial_id'])
    if x[0][0] != selected_trial_number:
        raise ValueError


    # =========================================================================
    # Upper grid
    # =========================================================================
    ax = plt.axes([0.3 ,0.4, 0.4, 0.5])
    #subplot2grid((6, 3), (0, 1), 3, 1)

    cmap = plt.cm.get_cmap('hsv')

    image = wave_plots.plot_array_col_grad(
        ax,
        phase_frame_cut[frame_idx, :].magnitude,
        np.ma.array(
            cohphgr_frame_cut[frame_idx, :].magnitude,
            mask=(results_prop['good_el'] is False)) /
        np.ma.array(
            np.abs(cohphgr_frame_cut[frame_idx, :].magnitude) * 2.,
            mask=(results_prop['good_el'] is False)),
        results_prop['good_el'], cmap)

    cbar = plt.colorbar(image)
 
    x_txt = cbar.ax.get_xlim()[0] + \
        (cbar.ax.get_xlim()[1] - cbar.ax.get_xlim()[0]) / 2.
    y_txt = cbar.ax.get_ylim()[0] + \
        (cbar.ax.get_ylim()[1] - cbar.ax.get_ylim()[0]) / 2.
    cbar.ax.text(
        x_txt, y_txt, 'phase (rad)',
        size=myaxissize, color='k', ha='center', va='center', rotation=90.)
    cpos=cbar.ax.get_position().bounds
    cbar.ax.set_position([cpos[0], cpos[1]+.05, cpos[2],cpos[3]-.1])

    plt.title("Phase and norm. gradient coherence", size=mylabelsize)


    # =========================================================================
    # complete signal
    # =========================================================================
    ax = plt.subplot2grid((6, 1), (4, 0), 1, 1)
    plt.plot([t_i - tpre_ms, t_i - tpre_ms], [4.5, -4], 'k', linewidth=3)
    plt.plot(
        signal_cut.times.rescale(pq.ms).magnitude - tpre,
        signal_cut.magnitude, color='gray')
    plt.plot(
        signal_f_cut.times.rescale(pq.ms).magnitude - tpre,
        signal_f_cut.magnitude, color=[1, 0.2, 0.7])
    wave_plots.plot_precise_events(
        results_events['event_stat'][selected_monkey][
            selected_subsession_cnd]['all'],
        -6, 6,
        1, meanstyle='k:', stdstyle=None)
    eticks = []
    eticks_labels = []
    for ev in wave_plots.event_names:
        eticks.append(
            results_events['event_stat'][selected_monkey][
                selected_subsession_cnd]['all'][ev][0])
        eticks_labels.append(ev)
    ax.set_xticks(eticks)
    ax.set_xticklabels(eticks_labels, rotation=45, ha='right')
    ax.tick_params(labelsize=myaxissize, size=0)
    plt.xlim((-tpre, signal_f_cut.times[-1].rescale(pq.ms).magnitude - tpre))
    plt.ylim(-4, 6)
    plt.xlabel('time (ms)', size=mylabelsize)
    ax.tick_params(labelsize=myaxissize)
    ax.set_yticks([-4, -2, 0, 2, 4])
    plt.axis('off')
    plt.title("LFP and $\\beta$-filtered LFP (single el.)", size=mylabelsize)

    # =========================================================================
    # Sketch
    # =========================================================================
#     ax = plt.subplot2grid((3, 2), (2, 1), 1, 1)
#     cmap = matplotlib.colors.ListedColormap([
#         "black",
#         wave_plots.pattern_colors[wave_main.state_names[
#             int(th_clusters_cut.magnitude[frame_idx])]],
#         'white'])
#     pic = matplotlib.image.imread(
#         './material/cartoons/pattern_' +
#         wave_main.state_names[
#             int(th_clusters_cut.magnitude[frame_idx])] + '.eps')
#     ax.imshow(pic[:, :, 0], aspect='auto', cmap=cmap)
#     ax.set_aspect('equal')
#     plt.axis('off')

#     plt.figtext(
#         0.6, 0.4,
#         "Classification:\n%s" %
#         wave_main.state_names[int(th_clusters_cut.magnitude[frame_idx])],
#         color=wave_plots.pattern_colors[wave_main.state_names[
#             int(th_clusters_cut.magnitude[frame_idx])]],
#         fontsize=20)
# 
#     plt.axis('off')

#     ax = plt.axes([0.02, 0.15, .1, .1])
#     pic = matplotlib.image.imread('./material/wave_ms_fig3/ampl.png')
#     ax.imshow(pic, aspect='auto')
#     ax.set_aspect('equal')
#     plt.axis('off')

# =============================================================================
# Initialize project
# =============================================================================

# Fetch job ID
job_id = wave_main.get_jobid()

# Create job(s) and prefetch files if required
jobparameters = wave_main.create_file_filter_task(
    wave_main.selected_datasets, wave_main.selected_filter_names,
    fetch=False)

selected_subsession = jobparameters[0][0]
selected_filter = wave_main.selected_filter_names[0]
selected_trial_id = 45

# selected_subsession = wave_main.selected_datasets[45]
# selected_filter = wave_main.selected_filter_names[1]
# selected_trial_id = 45

selected_subsession_name = selected_subsession[0]
selected_subsession_cnd = selected_subsession[1]

selected_monkey = wave_main.monkey_names[selected_subsession[0][0]]

print("Locating, fetching and verifying data")
data_dir, metadata_dir = wave_main.get_data_dir(selected_subsession_name)

# =============================================================================
# Initialize project
# =============================================================================

# Result directory contains filename, filter and trial to overcome problems
# with number of files in a directory
pc = projctrl.ProjectControl(
    project_name='collab-waves',
    script_name='ms_figs/movs1_small/movie_%s_filter_%s_trialid_%i' %
    (selected_subsession_name, selected_filter, selected_trial_id),
    clear_data=True)

# =============================================================================
# Set and save parameters
# =============================================================================

param_prop = h5pyw.load_h5(
    pc.result_domain + '/calc_waveproperties/' + selected_subsession_name +
    '_filter_' + selected_filter + '_param.h5')

param = {}

param['trigger'] = 'TS'
param['tpre'] = 200. * pq.ms
param['tpost'] = 3400. * pq.ms
param['trial_id'] = selected_trial_id
param['step_size'] = 2 * pq.ms
param['min_stat_duration'] = 5

mylabelsize = 14
myaxissize = 10

# =============================================================================
# Load data
# =============================================================================

print("Loading data")
subsessionobj = ReachGraspIO(
    os.path.join(data_dir, selected_subsession_name),
    metadata_dir=metadata_dir, print_diagnostic=False)

# Check if trial ID is correct and get linear trial number in file
# TODO: Remove if cutting and selection routines are no longer required.
trialids = subsessionobj.get_trial_ids(performance_codes=127)
try:
    selected_trial_number = trialids.index(selected_trial_id)
except:
    raise ValueError(
        "Selected trial ID does not correspond to a correct trial.")
param['trial_nr'] = selected_trial_number

# Read in data and create z-score data
neo_block = subsessionobj.read_block(
    n_starts=[None], n_stops=[None],
    channel_list=range(1, 97), nsx=2, events=True)
pick.map_as(
    neo_block, elephant.signal_processing.zscore, 
    annotations=None, inplace=False)

# Create filtered and z-scored data
neo_block_f = subsessionobj.read_block(
    n_starts=[None], n_stops=[None],
    channel_list=range(1, 97), nsx=2, events=True)
pick.map_as(
    neo_block_f, elephant.signal_processing.butter,
    annotations=None,
    lowpass_freq=param_prop['highcut'],
    highpass_freq=param_prop['lowcut'],
    order=param_prop['order'])

pick.map_as(
    neo_block_f, elephant.signal_processing.zscore, 
    annotations=None, inplace=False)

# Add trial epochs to the data
subsessionobj.add_trials_around(
    neo_block, name='mytrial', aligned_trigger=param['trigger'],
    tpre=param['tpre'], tpost=param['tpost'],
    aligned_offset=0 * pq.ms)
subsessionobj.add_trials_around(
    neo_block_f, name='mytrial', aligned_trigger=param['trigger'],
    tpre=param['tpre'], tpost=param['tpost'],
    aligned_offset=0 * pq.ms)

# Read neo frames
results_prop = h5pyw.load_h5(
    pc.result_domain + os.path.sep + 'calc_waveproperties' + os.path.sep +
    selected_subsession_name + '_filter_' + selected_filter + '_info.h5')

ho = nh5.NeoHdf5IO(
    filename=pc.result_domain + os.path.sep + 'calc_waveproperties' +
    os.path.sep + selected_subsession_name + '_filter_' + selected_filter +
    '_neo_frames.h5')
print("Loading phase frames")
phase_frames = ho.get(results_prop['phase_frames_path'], cascade=False)
print("Loading amplitude frames")
amp_frames = ho.get(results_prop['amp_frames_path'], cascade=False)
print("Loading phase gradient frames")
phgr_frames = ho.get(results_prop['phgr_frames_path'], cascade=False)
print("Loading coherence gradients")
cohphgr_frames = ho.get(results_prop['cohphgr_frames_path'], cascade=False)
ho.close()


# Choose an example electrode for the plot (i.e., first good linear electrode)
# param['sample_el'] = np.where(results_prop['good_el'])[0][0] + 1
param['sample_el'] = 50
if not results_prop['good_el'][param['sample_el']-1]:
    raise ValueError('Sample Electrode is not good')

# Cut signals to selected trial
# TODO: Unify into a single cutting routine
# TODO: Cutting seems inconsistent with respect to where zero is!
print("Cutting")

# The following return a signal of a specific trial cut to a specified trial
# ID, and times starting at 0.
signal_cut = pick.get_trials_AnalogSignal(
    neo_block,
    {'trial_id': [param['trial_id']]}, {'ca_id': [param['sample_el']]},
    reset_times=True)[0][1][0]

signal_f_cut = pick.get_trials_AnalogSignal(
    neo_block_f,
    {'trial_id': [param['trial_id']]}, {'ca_id': [param['sample_el']]},
    reset_times=True)[0][1][0]

# The following return a cut AnalogSignalArray where the i-th trial in the
# epocharray is cut out, and times start at 0.
phase_frame_cut = wave_main.cut_asa(
    neo_block.segments[0].epocharrays[0],
    phase_frames, reset_times=True, i=selected_trial_number)

amp_frame_cut = wave_main.cut_asa(
    neo_block.segments[0].epocharrays[0],
    amp_frames, reset_times=True, i=selected_trial_number)

phgr_frame_cut = wave_main.cut_asa(
    neo_block.segments[0].epocharrays[0],
    phgr_frames, reset_times=True, i=selected_trial_number)

cohphgr_frame_cut = wave_main.cut_asa(
    neo_block.segments[0].epocharrays[0],
    cohphgr_frames, reset_times=True, i=selected_trial_number)

# Load pooling parameters and results of power analysis
results_events = h5pyw.load_h5(
    pc.result_domain + '/calc_events/' + 'events.h5')


# =============================================================================
# Plot trial
# =============================================================================

print("Plotting")

# Base file name
filename = '%s_filter_%s_trialid_%i' % \
    (selected_subsession_name, selected_filter, selected_trial_id)

# Make individual frames using the frame maker
wave_plots.create_frames(
    frame_func=phasemovie_frames,
    num_frames=int(
        ((signal_cut.t_stop - signal_cut.t_start) /
            param['step_size']).simplified.magnitude),
    frame_interval=param['step_size'],
    pc=pc,
    frame_name=filename,
    figsize=(10.24, 7.80),
    keep_eps=False,
    spare_core=-1)

# Assemble movie
# wave_plots.movie_maker(
#     frames_directory=pc.figure_path + os.path.sep + "png",
#     movie_directory=pc.figure_path + os.path.sep + "mov",
#     frame_format="png",
#     movie_name=filename,
#     frames_per_sec=15,
#     quality=23,
#     scale_x=1024,
#     scale_y=780)

# Write parameters to disk
print("Saving  parameters - %s" % selected_subsession_name)
h5pyw.add_to_h5(
    pc.result_path + selected_subsession_name +
    '_filter_' + selected_filter + '_trialid_' +
    str(selected_trial_id) + '_param.h5', param,
    write_mode='w', overwrite_dataset=True)
