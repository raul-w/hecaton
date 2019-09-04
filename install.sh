#!/bin/bash

set -e

# set up conda environments
conda env create -f docker/environment_py3.yml && \
conda env create -f docker/environment_py2.yml && \

# install all other dependencies
mkdir hecaton_deps && \
cd hecaton_deps && \
wget https://github.com/PapenfussLab/gridss/releases/download/v2.0.1/gridss-2.0.1-gridss-jar-with-dependencies.jar && \
echo "export GRIDSS_JAR=$PWD/gridss-2.0.1-gridss-jar-with-dependencies.jar" >> ~/.bashrc && \
wget https://github.com/broadinstitute/picard/releases/download/2.18.23/picard.jar && \
echo "export PICARD=$PWD/picard.jar" >> ~/.bashrc && \
source activate hecaton_py2 && \
git clone --recursive https://github.com/hall-lab/speedseq && \
cd speedseq && \
make align && \
make sv && \
make config && \
echo "export PATH=$PATH:$PWD/bin" >> ~/.bashrc && \
source deactivate && \
cd ../.. && \
source ~/.bashrc
