#!/bin/bash
# Install code repo
git clone https://1662a097532b110aba67551f7f13273355c52efe@github.com/INM-6/reachgrasp-spikewave.git

# Create environment
virtualenv rgsw-env
source activate rgsw-env
pip install -r reachgrasp-spikewave/src/requirements.txt
