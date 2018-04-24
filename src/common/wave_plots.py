"""
Module containing global plotting utility functions for the wave project.
"""

import os.path
import subprocess
import multiprocessing

import numpy as np
import matplotlib.pyplot as plt

import quantities as pq

import wave_main as wave_main

# Colors of pattern types
cmap = [
    'DarkGray', 'Navy', 'DarkTurquoise',
    'YellowGreen', 'Orange', 'Crimson']

# Alternative
# (unclassified is black, then add palette for patterns)
# cmap = [(0, 0, 0, 1)]
# cmap.extend([plt.cm.get_cmap('gist_ncar', 20)(i) for i in [
#     2, 8, 12, 14, 16]])

# cmap.extend([plt.cm.get_cmap('gist_ncar', 16)(i) for i in range(1, 16, 3)])
# cmap.extend([plt.cm.get_cmap('nipy_spectral', 7)(i) for i in range(1, 6)])
# cmap.extend([plt.cm.get_cmap('hsv', 7)(i) for i in range(1, 6)])

pattern_colors = dict(zip(
    wave_main.state_names,
    cmap))

# Colors of measurements
meas_colors = {
    'meas_var_phase': 'Black',
    'meas_var_norm_phgr': 'Black',
    'meas_avg_cohphgr': 'Black',
    'meas_contingency_cohphgr': 'Black',
    'meas_circlishness_cohphgr': 'Black',
    'meas_outwardishness_cohphgr': 'Black'}

# Alternative
# meas_colors = {
#     'meas_var_phase': 'Indigo',
#     'meas_var_norm_phgr': 'MediumSlateBlue',
#     'meas_avg_cohphgr': 'SeaGreen',
#     'meas_contingency_cohphgr': 'YellowGreen',
#     'meas_circlishness_cohphgr': 'Orange',
#     'meas_outwardishness_cohphgr': 'Crimson'}

monkey_names_plot = {
    'Nikos2': 'monkey N',
    'Lilou': 'monkey L',
    'Tanya': 'monkey T'}


event_names = [
    'TS',
    'WS',
    'CUE-ON',
    'CUE-OFF',
    'GO',
    'SR',
    'RW']


def plot_precise_events(
        event_struct, ymin, ymax, scale,
        meanstyle='g-', stdstyle='g--', minmaxstyle=None):
    for en in event_names:
        m = event_struct[en][0] / scale
        s = event_struct[en][1] / scale
        m1 = event_struct[en][2] / scale
        m2 = event_struct[en][3] / scale

        if meanstyle:
            plt.plot([m, m], [ymin, ymax], meanstyle)
        if stdstyle:
            plt.plot([m - s, m - s], [ymin, ymax], stdstyle)
            plt.plot([m + s, m + s], [ymin, ymax], stdstyle)
        if minmaxstyle:
            plt.plot([m1, m1], [ymin, ymax], minmaxstyle)
            plt.plot([m2, m2], [ymin, ymax], minmaxstyle)


# Deprecated!
def plot_events(frame_interval, ymin, ymax):
    t = pq.Quantity(
        0.0, units=frame_interval.dimensionality).rescale(pq.ms).magnitude
    plt.plot([t, t], [ymin, ymax], 'g--')

    t = pq.Quantity(
        400.0, units=frame_interval.dimensionality).rescale(pq.ms).magnitude
    plt.plot([t, t], [ymin, ymax], 'g--')

    t = pq.Quantity(
        800.0, units=frame_interval.dimensionality).rescale(pq.ms).magnitude
    plt.plot([t, t], [ymin, ymax], 'g--')

    t = pq.Quantity(
        1100.0, units=frame_interval.dimensionality).rescale(pq.ms).magnitude
    plt.plot([t, t], [ymin, ymax], 'g--')

    t = pq.Quantity(
        2100.0, units=frame_interval.dimensionality).rescale(pq.ms).magnitude
    plt.plot([t, t], [ymin, ymax], 'g--')


def plot_array_col_grad(
        ax, colordata, gradientdata, el_mask, cmap, minc=-np.pi, maxc=np.pi,
        rot=0):
    """Plots a typical wave plot with colors and arrows

    Parameters
    ==========
    ax : int
        Handle to axis in which to plot the array.
    colordata : neo.AnalogSignalArray
        Contains the signal array of 100 signals to plot on the array.
        Convention: Signal 0 is the lower left, signal 9 is the lower right,
        signal 99 is the top right.
    gradientdata : neo.AnalogSignalArray
        Gradients to plot as arrow. Caution: Arrows point opposite to the
        gradient! If set to None, no arrows are plotted.
    el_mask : numpy.array of bool
        100x1 array of booleans that indicate whether to plot data from each of
        the 100 electrodes. True will plot the data, False will plot a black
        square.
    cmap : matplotlob.pyplot.cm.get_cmap
        Colormap to for plotting the phases.
    rot : float
        Rotates the array CW by rot degree. Currently partly supported.
        Default: 0.


    Returns
    =======
    numpy.array
        Image returned by pcolormesh.
    """
    # Broken electrodes
    cmap.set_bad(color='k')

    # Arrange data on a grid Note: -1 factor for phase gradients to show anti-
    # gradients that resemble the direction of propagation
    phase_grid = np.ma.array(np.reshape(
        colordata, (10, 10)), mask=(np.equal(el_mask, False)))

    # Plot phases as pcolormesh
    X, Y = np.meshgrid(np.arange(0, 11), np.arange(0, 11))
    if rot != 0:
        # rotate clockwise
        rot = 2.*np.pi-rot

        # affine transform
        X = X - 4.5
        Y = Y - 4.5
        Xr = X * np.cos(rot) - Y * np.sin(rot)
        Yr = X * np.sin(rot) + Y * np.cos(rot)
        X = Xr + 4.5
        Y = Yr + 4.5
    image = ax.pcolormesh(
        X, Y, phase_grid, vmin=minc, vmax=maxc, cmap=cmap)

    # Plot arrows
    if gradientdata is not None:
        direc_grid = -np.ma.array(np.reshape(
            gradientdata, (10, 10)), mask=(el_mask is False))
        largedirec_grid = np.ma.array([
            [np.mean(direc_grid[0:5, 0:5]),
                np.mean(direc_grid[0:5, 5:10])],
            [np.mean(direc_grid[5:10, 0:5]),
                np.mean(direc_grid[5:10, 5:10])]]) * 5.0

        # Scale factor 10 to account for 10x10 grid
        X, Y = np.meshgrid(np.arange(0.5, 10.5), np.arange(0.5, 10.5))
        ax.quiver(
            X, Y, np.real(direc_grid), np.imag(direc_grid),
            color='k', units='x', scale=1.0, scale_units='x', pivot='middle')

        # Plot big arrows
        X, Y = np.meshgrid([2.5, 7.5], [2.5, 7.5])
        ax.quiver(
            X, Y, np.real(largedirec_grid), np.imag(largedirec_grid),
            color='w', units='x', scale=1.0, scale_units='x', pivot='middle',
            width=0.2, headaxislength=5)

    # Adapt axis
    ax.set_aspect('equal', 'box-forced', 'C')
    ax.tick_params(
        left='off', bottom='off', top='off', right='off', labelsize='xx-small')
    ax.set_xticks(np.arange(0.5, 10.5, 1))
    ax.set_xticklabels(
        ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
    ax.set_yticks(np.arange(0.5, 10.5, 1))
    ax.set_yticklabels(
        ['1', '11', '21', '31', '41', '51', '61', '71', '81', '91'])

    return image


def _frame_worker(
        frame_func, frame_i, t_i, frame_interval, num_frames,
        frame_name, pc, figsize, keep_eps):
    """
    Helper function for create_frames to perform the plotting and saving of
    each frame
    """

    print(
        "Creating Frame %i/%i at %.1f ms" %
        (frame_i, num_frames, t_i.rescale(pq.ms).magnitude))

    # Create new figure for frame
    if figsize is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=figsize)

    # Draw frame
    frame_func(frame_i, t_i, frame_interval)

    # Add frame number to filename and save figure
    fname = frame_name + '_' + str(frame_i).zfill(5)
    pc.save_fig(
        fname,
        save_png=True, save_jpg=False, save_eps=keep_eps, save_pdf=False)
    plt.close(fig)


def create_frames(
        frame_func, num_frames, frame_interval,
        pc, frame_name, figsize=None, keep_eps=False, spare_core=-1):
    """
    Creates and saves all the frames of a movie sequence.

    Parameters
    ==========
    frame_func : function
        Function that draws a single frame by filling the current figure. The
        function must be of the form frame_func(frame_i, t_i, frame_interval),
        where the function plots the frame number frame_i, time t_i (starting
        at zero, in units of frame_interval), at a frame interval of
        frame_interval.
    num_frames : int
        Number of frames to plot.
    frame_interval : Quantity
        Time interval between two consecutive movie frames.
    pc : projctrl.ProjectControl
        Project control structure used to save each movie frame.
    frame_name : str
        Base name of each movie frame.
    figsize : tuple
        Passed to figure. If None is specified, the default is used.
    keep_eps : bool
        If true, keeps the eps figure.
    spare_core : int
        Number of cores to spare in multiprocessing. For example, if set to 2,
        only 6 cores on a 8 core machines will not be used. If -1, no
        multiprocessing is used at all. Default: -1

    """

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
        for frame_i in range(num_frames):
            # Time index of frame
            t_i = frame_i * frame_interval
            _frame_worker(
                frame_func, frame_i, t_i, frame_interval, num_frames,
                frame_name, pc, figsize, keep_eps)
    else:
        pool = multiprocessing.Pool(processes=num_core)

        # Launch jobs
        results = []
        for frame_i in range(num_frames):
            # Time index of frame
            t_i = frame_i * frame_interval
            print(
                "Submitting Frame %i/%i at %.1f ms" %
                (frame_i, num_frames, t_i.rescale(pq.ms).magnitude))
            results.append(pool.apply_async(
                func=_frame_worker,
                args=[
                    frame_func, frame_i, t_i, frame_interval, num_frames,
                    frame_name, pc, figsize, keep_eps]))
        pool.close()
        for r in results:
            r.get()
        print("Multiprocessing done.")
        pool.join()
        print("Multiprocessing joined.")
        pool.terminate()


def movie_maker(frames_directory, movie_directory, frame_format, movie_name,
                frames_per_sec=20, quality=23, scale_x=1920, scale_y=1080):
    """
    Makes a movie from a given frame series.

    Parameters
    ==========
    frames_directory : str
        Directory where all frames of the movie are located
    movie_directory : str
        Directory in which movie should be saved
    frame_format : str
        Frame format as string (e.g. 'png')
    movie_name : str
        Filename of the movie (without extension)
    frames_per_sec : int
        Number of frames per second. Default: 20
    quality
        This is the crf argument of avconv. Here, higher is smaller file size
        (0=lossless, 23=average, 51=worst quality). Default: 23
    scale_x, scale_y : int
        X and Y resolution of movie. Default: Full HD (1920x1080)
    """

    # mp4 using avconv
    command = (
        'avconv',
        '-y',  # Overwrite existing file
        '-v', 'quiet',
        '-r', str(frames_per_sec),  # Frame rate
        '-i', frames_directory + os.path.sep + \
        movie_name + '_%05d.' + frame_format,
        '-c:v', 'h264',  # MP4 Codec
        '-crf', str(quality),  # Quality
        '-vf', 'scale=' + str(scale_x) + ':' + str(scale_y),  # Resolution
        movie_directory + os.path.sep + movie_name + '.mp4')

    print("\n\nExecuting:\n %s\n\n" % ' '.join(command))
    subprocess.check_call(command)
