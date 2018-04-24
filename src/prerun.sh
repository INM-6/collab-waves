#!/bin/bash
# Install code repo
git clone https://git@github.com/INM-6/reachgrasp-spikewave.git

# Create environment
virtualenv rgsw-env
source activate rgsw-env
pip install -r reachgrasp-spikewave/src/requirements.txt
