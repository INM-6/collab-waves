import glob
import os.path as path
import os

import sys
sys.path.insert(0,"src")
sys.path.insert(0,"src/common")

from wave_main import selected_datasets

datapath = {
    'l': "../DataLilou", 
    'i': "../DataNikos2"}

rule generate_events:
    input:
        expand(
            "datasets/{session}.nev",
            session=[
                x[0] for x in selected_datasets])
    output:
        "ProjectsData/reachgrasp-spikewave/results/calc_trial_events/events.h5"
    shell:
        "python reachgrasp-spikewave/src/calc_trial_events ;"


rule generate_phases:
    input:
        lambda wildcards: "datasets/" + wildcards.session + ".nev"
    output:
        "ProjectsData/reachgrasp-spikewave/results/calc_waveproperties/{session}_filter_beta_large_neo_frames.h5"
    run:
        i = None
        for j,k in enumerate(selected_datasets):
            if j[0] == wildcards.session:
                i = k
        if i==None:
            raise ValueError("Unknown file.")
	    shell("python reachgrasp-spikewave/src/calc_waveproperties " + str(i))

rule generate_bollywood:
    input:
        lambda wildcards: "ProjectsData/reachgrasp-spikewave/results/calc_waveproperties/"+ wildcards.session + "_filter_beta_large_neo_frames.h5"
    output:
        "ProjectsData/reachgrasp-spikewave/figs/ms_figs/movs1_small/movie_{session}_filter_beta_large_trialid_45/png/{session}_filter_beta_large_trialid_45_00000.png"
    run:
        i = None
        for j,k in enumerate(wave_main.selected_datasets):
            if j[0] == wildcards.session:
                i = k
        if i != None:
	        shell("python reachgrasp-spikewave/src/calc_waveproperties " + str(i))
        else:
            print("Unknown file.")


rule all:
    input:
        expand(
	    "ProjectsData/reachgrasp-spikewave/figs/ms_figs/movs1_small/movie_{session}_filter_beta_large_trialid_45/png/${session}_filter_beta_large_trialid_45_00000.png",
            session=[ 
                x[0] for x in selected_datasets])
