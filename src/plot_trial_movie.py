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

import jelephant.core.pick
from rg.rgio import ReachGraspIO

import projctrl
import h5py_wrapper.wrapper

import common.wave_main as wave_main
import common.wave_plots as wave_plots

__updated__ = "2016-05-22"


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
        print x[0][0]
        print selected_trial_number
        raise ValueError

    plt.figtext(
        .5, 0.96,
        'Session: %s, Filter: %.1f-%.1f Hz, Trial ID: %i, Type: %s, '
        'Frame %i [%.1f ms]' % (
            selected_subsession_name,
            param_prop['lowcut'], param_prop['highcut'],
            selected_trial_id,
            subsessionobj.trial_type_str[
                neo_block.segments[0].epocharrays[0].annotations[
                    'trial_type'][selected_trial_number]],
            frame_i, t_i - tpre_ms),
        size=24., ha='center', va='center')

    # =========================================================================
    # Upper grid
    # =========================================================================
    ax = plt.subplot2grid((14, 3), (0, 0), 5, 1)

    cmap = plt.cm.get_cmap('hsv')

    image = wave_plots.plot_array_col_grad(
        ax,
        phase_frame_cut[frame_idx, :].magnitude,
        np.ma.array(
            phgr_frame_cut[frame_idx, :].magnitude,
            mask=(results_prop['good_el'] is False)) /
        np.ma.array(
            np.abs(phgr_frame_cut[frame_idx, :].magnitude) * 2.,
            mask=(results_prop['good_el'] is False)),
        results_prop['good_el'], cmap)

    cbar = plt.colorbar(image)
    x_txt = cbar.ax.get_xlim()[0] + \
        (cbar.ax.get_xlim()[1] - cbar.ax.get_xlim()[0]) / 2.
    y_txt = cbar.ax.get_ylim()[0] + \
        (cbar.ax.get_ylim()[1] - cbar.ax.get_ylim()[0]) / 2.
    cbar.ax.text(
        x_txt, y_txt, 'phase (rad)',
        size=6., color='k', ha='center', va='center', rotation=90.)

    plt.title("Phases and norm. phase gradients ")

    # =========================================================================
    # lower grid
    # =========================================================================
#     ax = plt.subplot2grid((14, 3), (6, 0), 5, 1)
#
#     cma = plt.cm.get_cmap('hsv')(range(256))
#     for i in range(cma.shape[0]):
#         cmh = np.array(colorsys.rgb_to_hsv(cma[i, 0], cma[i, 1], cma[i, 2]))
#         cmh[1] = meas['meas_avg_amp'][frame_idx] / \
#             np.max(meas['meas_avg_amp']) * 2.0 / 3.0 + 0.32
#         cma[i, 0], cma[i, 1], cma[i, 2], = \
#             colorsys.hsv_to_rgb(cmh[0], cmh[1], cmh[2])
#     cmap = matplotlib.colors.ListedColormap(cma)
#
#     image = wave_plots.plot_array_col_grad(
#         ax,
#         phase_frame_cut[frame_idx, :].magnitude,
#         cohphgr_frame_cut[frame_idx, :].magnitude,
#         results_prop['good_el'], cmap)
#
#     cbar = plt.colorbar(image)
#     x_txt = cbar.ax.get_xlim()[0] + \
#         (cbar.ax.get_xlim()[1] - cbar.ax.get_xlim()[0]) / 2.
#     y_txt = cbar.ax.get_ylim()[0] + \
#         (cbar.ax.get_ylim()[1] - cbar.ax.get_ylim()[0]) / 2.
#     cbar.ax.text(
#         x_txt, y_txt, 'phase (rad)',
#         size=6., color='k', ha='center', va='center', rotation=90.)
#
#     plt.title("Phases, ampl. and coherence phase gradients")

    ax = plt.subplot2grid((14, 3), (6, 0), 5, 1)

    cmap = plt.cm.get_cmap('jet')

    image = wave_plots.plot_array_col_grad(
        ax,
        amp_frame_cut[frame_idx, :].magnitude,
        None,
        results_prop['good_el'], cmap, minc=0, maxc=4)

    cbar = plt.colorbar(image)
    x_txt = cbar.ax.get_xlim()[0] + \
        (cbar.ax.get_xlim()[1] - cbar.ax.get_xlim()[0]) / 2.
    y_txt = cbar.ax.get_ylim()[0] + \
        (cbar.ax.get_ylim()[1] - cbar.ax.get_ylim()[0]) / 2.
    cbar.ax.text(
        x_txt, y_txt, 'amplitude',
        size=6., color='k', ha='center', va='center', rotation=90.)

    plt.title("Amplitudes")

    # =========================================================================
    # complete signal
    # =========================================================================
    ax = plt.subplot2grid((14, 3), (12, 0), 2, 1)
    for clus in range(5):
        for z in range(len(th_classepochs[clus].times)):
            t1 = th_classepochs[clus].times[z].rescale(pq.ms).magnitude
            t2 = t1 + \
                th_classepochs[clus].durations[z].rescale(pq.ms).magnitude
            plt.fill_between(
                [t1 - tpre, t2 - tpre], -4, 4,
                color=wave_plots.pattern_colors[wave_main.state_names[clus + 1]], alpha=1., lw=0)
    plt.plot(
        signal_cut.times.rescale(pq.ms).magnitude - tpre,
        signal_cut.magnitude, 'k-')
    plt.plot(
        signal_f_cut.times.rescale(pq.ms).magnitude - tpre,
        signal_f_cut.magnitude, 'r-')
    plt.plot([t_i - tpre_ms, t_i - tpre_ms], [4, -4], 'r')
    wave_plots.plot_events(frame_interval, -4, 4)
    plt.xlim((-tpre, signal_f_cut.times[-1].rescale(pq.ms).magnitude - tpre))
    plt.ylim(-4, 4)
    plt.xlabel('t (ms)')
    ax.tick_params(labelsize='xx-small')

    # =========================================================================
    # rolling signal
    # =========================================================================
    ax = plt.subplot2grid((15, 3), (0, 1), 2, 1)

    for clus in range(5):
        for z in range(len(th_classepochs[clus].times)):
            t1 = th_classepochs[clus].times[z].rescale(pq.ms).magnitude
            t2 = t1 + \
                th_classepochs[clus].durations[z].rescale(pq.ms).magnitude
            plt.fill_between(
                [t1 - tpre, t2 - tpre], -4, 4,
                color=wave_plots.pattern_colors[wave_main.state_names[clus + 1]], alpha=1., lw=0)
    plt.plot(
        signal_cut.times.rescale(pq.ms).magnitude - tpre,
        signal_cut.magnitude, 'k-')
    plt.plot(
        signal_f_cut.times.rescale(pq.ms).magnitude - tpre,
        signal_f_cut.magnitude, 'r-')
    plt.plot(
        phase_frame_cut.times.rescale(pq.ms).magnitude - tpre,
        phase_frame_cut[:, param['sample_el'] - 1].magnitude / 2,
        '--', color=[0.5, 0.5, 0.5])
    plt.plot(
        amp_frame_cut.times.rescale(pq.ms).magnitude - tpre,
        amp_frame_cut[:, param['sample_el'] - 1].magnitude,
        '--', color=[0, 0.5, 0])
    plt.plot(
        amp_frame_cut.times.rescale(pq.ms).magnitude - tpre,
        - amp_frame_cut[:, param['sample_el'] - 1].magnitude,
        '--', color=[0, 0.5, 0])
    plt.plot(
        [t_i - tpre_ms, t_i - tpre_ms], [4, -4], 'r')
    wave_plots.plot_events(frame_interval, -4, 4)
    plt.title('Signal, phase and amplitude')
    plt.xlim((
        (t_i - 400 * pq.ms).rescale(pq.ms).magnitude - tpre,
        (t_i + 400 * pq.ms).rescale(pq.ms).magnitude - tpre))
    plt.ylim(-4, 4)
    plt.xlabel('t (ms)')
    ax.tick_params(labelsize='xx-small')

    # =========================================================================
    # classification centroids
    # =========================================================================
#     ax = plt.subplot2grid((15, 3), (3, 1), 3, 1)
#     cmap = plt.cm.get_cmap('jet', len(meas))
#
#     # only print first 7 measures used in classification
#     for i in range(7):
#         plt.bar(np.arange(5) + i * 0.1, km_centroids[:, i], width=0.1,
#         color=cmap(i), edgecolor=None)
#     plt.xlim(0, 5)
#     plt.xlabel('Detected State')
#     ax.set_xticks(np.arange(0.4, 5.4, 1))
#     ax.set_xticklabels(['1', '2', '3', '4', '5'])
#     ax.tick_params(labelsize='xx-small')

    # =========================================================================
    # PCA
    # =========================================================================
#     ax = plt.subplot2grid((15, 3), (9, 1), 6, 1)
#     cmap = plt.cm.get_cmap('jet', 5)
#     plt.plot(pca_meas[:frame_idx, 0], pca_meas[:frame_idx, 1],
#         'o', c=[0.8, 0.8, 0.8], markersize=3)
#     plt.plot(pca_meas[frame_idx, 0], pca_meas[frame_idx, 1],
#         'ok', markersize=3)
#     for i in range(5):
#         # plot larger circle of the current state
#         if int(km_clusters_cut.magnitude[frame_idx]) == i:
#             ms = 25
#         else:
#             ms = 10
#         plt.plot(pca_centroid[i, 0], pca_centroid[i, 1],
#             'o', c=cmap(i), markersize=ms)
#     plt.xlim(-4, 4)
#     plt.ylim(-4, 4)
#     plt.xlabel('PC1')
#     plt.ylabel('PC2')
#     plt.title('k-means analysis')

    # =========================================================================
    # analysis-type measures
    # =========================================================================
    meas_str_plt = [
        'meas_avg_amp',
        'meas_avg_vel',
        'meas_avg_vel_dir',
        'meas_avg_dir']
    meas_name = [
        'Average amplitude',
        'Average phase velocity',
        'Average phase velocity in direction of gradient',
        'Average direction']
    cmap = plt.cm.get_cmap('jet', len(meas))
    for i, s in enumerate(meas_str_plt):
        ax = plt.subplot2grid((15, 3), (i * 2 + 4, 1), 1, 1)
        plt.plot(
            meas_cut[s].times.rescale(pq.ms).magnitude - tpre,
            meas_cut[s].magnitude,
            color=cmap(1.0 * (i + 10) / (len(meas_str) - 1)))
        plt.plot(
            [t_i - tpre_ms, t_i - tpre_ms],
            [np.min(meas_cut[s].magnitude), np.max(meas_cut[s].magnitude)],
            'r')
        wave_plots.plot_events(
            frame_interval, np.min(meas_cut[s].magnitude),
            np.max(meas_cut[s].magnitude))
        plt.title(meas_name[i])
        plt.xlim(
            (t_i - 400 * pq.ms).rescale(pq.ms).magnitude - tpre,
            (t_i + 400 * pq.ms).rescale(pq.ms).magnitude - tpre)
        plt.ylim(
            np.min(meas_cut[s].magnitude),
            np.max(meas_cut[s].magnitude))
        ax.tick_params(labelsize='xx-small')

    # =========================================================================
    # classification-type measures
    # =========================================================================
    meas_str_plt = [
        'meas_var_phase',
        'meas_var_norm_phgr',
        'meas_avg_cohphgr',
        'meas_contingency_cohphgr',
        'meas_circlishness_cohphgr',
        'meas_outwardishness_cohphgr']
    meas_name = [
        'Circular variance of phase $\sigma_p(t)$',
        'Circular variance of phase gradient $\sigma_g(t)$',
        'Average length of coh. phase grad. $\mu_c(t)$',
        'Average continuity $C(t)$',
        'Radial-perpendicular $R_{\perp}(t)$',
        'Radial-parallel $R_{\parallel}(t)$']
    #      #4B0082 #7B68EE #2E8B57 #9ACD32 #FFA500 #DC143C
    #     \newrgbcolor{meas1color}{0.2941 0.0    0.5098}
    #     \newrgbcolor{meas2color}{0.4842 0.4078 0.9333}
    #     \newrgbcolor{meas3color}{0.1804 0.5451 0.3412}
    #     \newrgbcolor{meas4color}{0.6040 0.8039 0.1961}
    #     \newrgbcolor{meas5color}{1.0000 0.6471 0.0000}
    #     \newrgbcolor{meas6color}{0.8627 0.0783 0.2353}
    meas_colors = [
        'Indigo', 'MediumSlateBlue', 'SeaGreen', 'YellowGreen',
        'Orange', 'Crimson']

    cmap = plt.cm.get_cmap('jet', len(meas))
    threshnr = 1

    for i, s in enumerate(meas_str_plt):
        ax = plt.subplot2grid(
            (len(meas_str_plt) * 2 - 1, 3),
            (i * 2, 2), 1, 1)

        # Plot axis
        plt.plot(
            meas_cut[s].times.rescale(pq.ms).magnitude - tpre,
            meas_cut[s].magnitude, color=meas_colors[i])

        # Plot marker of current position
        plt.plot([t_i - tpre_ms, t_i - tpre_ms], [-1, 1], 'r')

        # Plot marker for all events
        wave_plots.plot_events(pq.ms, -1, 1)

        # Draw thresholds
        for z in param_class['thresh'][i]:
            plt.axhline(
                z, meas_cut[s].times.rescale(pq.ms).magnitude[0] - tpre,
                meas_cut[s].times.rescale(pq.ms).magnitude[-1] - tpre,
                color='gray', linestyle='--')
            plt.text(
                (t_i + 410 * pq.ms).rescale(pq.ms).magnitude - tpre,
                z, r'$\theta_{' + str(threshnr) + '}$', color='gray')
            threshnr += 1
        plt.title(meas_name[i], color=meas_colors[i])
        plt.xlim(
            (t_i - 400 * pq.ms).rescale(pq.ms).magnitude - tpre,
            (t_i + 400 * pq.ms).rescale(pq.ms).magnitude - tpre)
        if i < 4:
            plt.ylim(0, 1)
        else:
            plt.ylim(-1, 1)
        ax.tick_params(labelsize='xx-small')

#     cmap = plt.cm.get_cmap('jet', len(meas))
#     for i, s in enumerate([meas_str[_]
#             for _ in [0, 1, 2, 3, 4, 5, 6, 7, 9, 11]]):
#         ax = plt.subplot2grid((19, 3), (i * 2, 2), 1, 1)
#         plt.plot(meas_cut[s].times.rescale(pq.ms).magnitude,
#             meas_cut[s].magnitude, color=cmap(1.0 * i / (len(meas_str) - 1)))
#         plt.plot([t_i, t_i],
#             [np.min(meas_cut[s].magnitude), np.max(meas_cut[s].magnitude)],
#             'r')
#         plot_events(frame_interval, np.min(meas_cut[s].magnitude),
#             np.max(meas_cut[s].magnitude))
#         plt.title("Property: " + s)
#         plt.xlim((t_i - 400 * pq.ms).rescale(pq.ms).magnitude,
#             (t_i + 400 * pq.ms).rescale(pq.ms).magnitude)
#         plt.ylim(np.min(meas_cut[s].magnitude),
#             np.max(meas_cut[s].magnitude))
#         ax.tick_params(labelsize='xx-small')

    # =========================================================================
    # Sketch
    # =========================================================================
    ax = plt.subplot2grid((15, 3), (12, 1), 3, 1)
    plt.title(
        "Threshold detection: %s" %
        wave_main.state_names[int(th_clusters_cut.magnitude[frame_idx])],
        color=wave_plots.pattern_colors[wave_main.state_names[
            int(th_clusters_cut.magnitude[frame_idx])]])
    cmap = matplotlib.colors.ListedColormap([
        "black",
        wave_plots.pattern_colors[wave_main.state_names[
            int(th_clusters_cut.magnitude[frame_idx])]],
        'white'])
    pic = matplotlib.image.imread(
        './material/cartoons/pattern_' +
        wave_main.state_names[
            int(th_clusters_cut.magnitude[frame_idx])] + '.eps')
    ax.imshow(pic[:, :, 0], aspect='auto', cmap=cmap)
    ax.set_aspect('equal')
    plt.axis('off')

    # TODO: Remove old picture code if new one works well
#     xy = (0, 0)
#     imagebox = matplotlib.offsetbox.OffsetImage(pic, zoom=.14)
#     ab = matplotlib.offsetbox.AnnotationBbox(
#         imagebox,
#         xy, xybox=(0., 0.), xycoords='data', pad=0.,
#         boxcoords="offset points", box_alignment=(0.5, 0.1))
#     ax.add_artist(ab)
    plt.axis('off')


# =============================================================================
# Initialize project
# =============================================================================

# Fetch job ID
job_id = wave_main.get_jobid()

# Create job(s) and prefetch files if required
jobparameters = wave_main.create_file_filter_task(
    wave_main.selected_datasets, wave_main.selected_filter_names,
    fetch=False)

selected_subsession = jobparameters[1]
selected_filter = jobparameters[1]

# Manual setting
selected_subsession = wave_main.selected_datasets[60]
selected_filter = wave_main.selected_filter_names[0]
selected_trial_id = 45

selected_subsession = wave_main.selected_datasets[1]
selected_filter = wave_main.selected_filter_names[0]
selected_trial_id = 81

selected_subsession = wave_main.selected_datasets[10]
selected_filter = wave_main.selected_filter_names[0]
selected_trial_id = 51

# selected_subsession = wave_main.selected_datasets[1]
# selected_filter = wave_main.selected_filter_names[1]
# selected_trial_id = 81

# selected_subsession = wave_main.selected_datasets[15]
# selected_filter = wave_main.selected_filter_names[0]
# selected_trial_id = 45

# selected_subsession = wave_main.selected_datasets[30]
# selected_filter = wave_main.selected_filter_names[1]
# selected_trial_id = 45

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
    project_name='reachgrasp-waves',
    script_name='plot_trial_movie/movie_%s_filter_%s_trialid_%i' %
    (selected_subsession_name, selected_filter, selected_trial_id),
    clear_data=True)

# =============================================================================
# Set and save parameters
# =============================================================================

param_prop = h5py_wrapper.wrapper.load_h5(
    pc.result_domain + '/calc_waveproperties/' + selected_subsession_name +
    '_filter_' + selected_filter + '_param.h5')
param_meas = h5py_wrapper.wrapper.load_h5(
    pc.result_domain + '/calc_wavemeas/' + selected_subsession_name +
    '_filter_' + selected_filter + '_param.h5')
param_class = h5py_wrapper.wrapper.load_h5(
    pc.result_domain + '/calc_waveclass/' + selected_subsession_name +
    '_filter_' + selected_filter + '_param.h5')

param = {}

param['trigger'] = 'TS'
param['tpre'] = 100. * pq.ms
param['tpost'] = 2700. * pq.ms
param['trial_id'] = selected_trial_id
param['step_size'] = 2 * pq.ms
param['min_stat_duration'] = 5


# =============================================================================
# Load data
# =============================================================================

print("Loading data")
subsessionobj = ReachGraspIO(
    os.path.join(data_dir, selected_subsession_name),
    metadata_dir=metadata_dir, print_diagnostic=False)

# Check if trial ID is correct and get linear trial number in file
# TODO: Remove if cutting and selection routines are no longer required.
trialids = subsessionobj.get_trial_ids(performance_codes=['correct_trial'])
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
jelephant.core.pick.map_as(
    neo_block, wave_main.applyzscore, annotations=None)

# Create filtered and z-scored data
neo_block_f = subsessionobj.read_block(
    n_starts=[None], n_stops=[None],
    channel_list=range(1, 97), nsx=2, events=True)
jelephant.core.pick.map_as(
    neo_block_f, wave_main.applyfilter,
    annotations=None, lowcut=param_prop['lowcut'],
    highcut=param_prop['highcut'], order=param_prop['order'])
jelephant.core.pick.map_as(
    neo_block_f, wave_main.applyzscore, annotations=None)

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
results_prop = h5py_wrapper.wrapper.load_h5(
    pc.result_domain + os.path.sep + 'calc_waveproperties' + os.path.sep +
    selected_subsession_name + '_filter_' + selected_filter + '_info.h5')
results_meas = h5py_wrapper.wrapper.load_h5(
    pc.result_domain + os.path.sep + 'calc_wavemeas' + os.path.sep +
    selected_subsession_name + '_filter_' + selected_filter + '_info.h5')
results_class = h5py_wrapper.wrapper.load_h5(
    pc.result_domain + os.path.sep + 'calc_waveclass' + os.path.sep +
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
# TODO: Remove or fix
# print("Loading phase latency frames")
# phlat_frames = ho.get(results_prop['phlat_frames_path'], cascade=False)
ho.close()

# TODO: Add in if required
# read neo measures
meas_str = [
    'meas_var_phase',
    'meas_var_phgr',
    'meas_var_norm_phgr',
    'meas_var_cohphgr',
    'meas_var_norm_cohphgr',
    'meas_avg_cohphgr',
    # 'meas_contingency_phgr',
    'meas_contingency_cohphgr',
    # 'meas_circlishness_phgr',
    'meas_circlishness_cohphgr',
    # 'meas_outwardishness_phgr',
    'meas_outwardishness_cohphgr',
    'meas_avg_amp',
    'meas_avg_vel',
    'meas_avg_vel_dir',
    'meas_avg_dir']
meas = {}
ho = nh5.NeoHdf5IO(
    filename=pc.result_domain + os.path.sep + 'calc_wavemeas' + os.path.sep +
    selected_subsession_name + '_filter_' + selected_filter +
    '_neo_measures.h5')
for i, s in enumerate(meas_str):
    print("Loading " + s)
    meas[meas_str[i]] = ho.get(
        results_meas[s + '_path'], cascade=False)
ho.close()

# Read Neo classifications
print("Loading classification")
ho = nh5.NeoHdf5IO(
    filename=pc.result_domain + os.path.sep + 'calc_waveclass' + os.path.sep +
    selected_subsession_name + '_filter_' + selected_filter + '_neo_class.h5')
th_clusters = ho.get(results_class['th_clusters_path'], cascade=False)
km_clusters = ho.get(results_class['km_clusters_path'], cascade=False)
hc_clusters = ho.get(results_class['hc_clusters_path'], cascade=False)
hc_linkageZ = ho.get(results_class['hc_linkage_path'], cascade=False)
km_centroids = results_class['km_centroids']
hc_lastdist = results_class['hc_lastdist']
ho.close()

# Choose an example electrode for the plot (i.e., first good linear electrode)
param['sample_el'] = np.where(results_prop['good_el'])[0][0] + 1

# Cut signals to selected trial
# TODO: Unify into a single cutting routine
# TODO: Cutting seems inconsistent with respect to where zero is!
print("Cutting")

# The following return a signal of a specific trial cut to a specified trial
# ID, and times starting at 0.
signal_cut = jelephant.core.pick.get_trials_AnalogSignal(
    neo_block,
    {'trial_id': [param['trial_id']]}, {'ca_id': [param['sample_el']]},
    reset_times=True)[0][1][0]

signal_f_cut = jelephant.core.pick.get_trials_AnalogSignal(
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

meas_cut = {}
for s in meas_str:
    meas_cut[s] = wave_main.cut_as(
        neo_block.segments[0].epocharrays[0],
        meas[s], reset_times=True, i=selected_trial_number)

th_clusters_cut = wave_main.cut_as(
    neo_block.segments[0].epocharrays[0],
    th_clusters, reset_times=True, i=selected_trial_number)

km_clusters_cut = wave_main.cut_as(
    neo_block.segments[0].epocharrays[0],
    km_clusters, reset_times=True, i=selected_trial_number)

th_classepochs, th_clusters_cut_reduced = \
    wave_main.create_class_epochs(
        th_clusters_cut, param['min_stat_duration'])

# TODO: Make sure +1 is not required
km_classepochs, km_clusters_cut_reduced = \
    wave_main.create_class_epochs(
        km_clusters_cut + 1, param['min_stat_duration'])

# construct PCA representation
# X = np.vstack([meas[i].magnitude for i in range(7)]).T
# whitened_X = vqcluster.whiten(X)
# X_cut = np.vstack([meas[i].magnitude for i in range(7)]).T
# whitened_X_cut = vqcluster.whiten(X_cut)
#
# n_components = 3
# pca = deco.PCA(n_components)
# pca.fit(whitened_X)
# pca_meas = pca.transform(whitened_X_cut)
# pca_centroid = pca.transform(km_centroids)


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
    figsize=(19.20, 10.80),
    keep_eps=False,
    spare_core=2)

# Assemble movie
wave_plots.movie_maker(
    frames_directory=pc.figure_path + os.path.sep + "png",
    movie_directory=pc.figure_path + os.path.sep + "mov",
    frame_format="png",
    movie_name=filename,
    frames_per_sec=15,
    quality=23,
    scale_x=1920,
    scale_y=1080)

# Write parameters to disk
print("Saving  parameters - %s" % selected_subsession_name)
h5py_wrapper.wrapper.add_to_h5(
    pc.result_path + selected_subsession_name +
    '_filter_' + selected_filter + '_trialid_' +
    str(selected_trial_id) + '_param.h5', param,
    write_mode='w', overwrite_dataset=True)
