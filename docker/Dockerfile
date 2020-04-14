FROM phusion/baseimage:0.11
MAINTAINER Raúl Wijfjes <raul.wijfjes@wur.nl>
LABEL authors="raul.wijfjes@wur.nl" \
    description="Docker image containing all dependencies of hecaton"

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Install basic applications
RUN apt-get update && \
	apt-get install -y curl gcc g++ wget git unzip tar gzip bzip2 xz-utils pigz make zlib1g-dev libncurses5-dev openjdk-8-jdk ant ca-certificates ca-certificates-java && \
	apt-get clean && \
	update-ca-certificates -f && \
	rm -rf /var/lib/apt/lists/* && \ 
	rm -rf /var/cache/oracle-jdk8-installer

# Install Miniconda (https://github.com/ContinuumIO/docker-images/blob/master/miniconda3/Dockerfile)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
	echo "conda activate base" >> ~/.bashrc

# Install Nextflow
RUN wget -qO- https://get.nextflow.io | bash && \
	chmod 755 nextflow && \
	mv nextflow /usr/bin/

# Install GRIDSS and Picard
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v2.0.1/gridss-2.0.1-gridss-jar-with-dependencies.jar && \
	wget https://github.com/broadinstitute/picard/releases/download/2.18.23/picard.jar

# Set path to GRIDSS jar file
ENV GRIDSS_JAR $PWD/gridss-2.0.1-gridss-jar-with-dependencies.jar

# Set path to Picard jar file
ENV PICARD $PWD/picard.jar

# Set path to Conda
ENV PATH /opt/conda/bin:$PATH 

# Install the rest of the required tools using Conda
COPY environment_py3.yml /
COPY environment_py2.yml /
RUN conda env create -f /environment_py3.yml && \
	conda env create -f /environment_py2.yml && \
	conda clean -a

# Install Speedseq
SHELL ["/bin/bash", "-c"]

RUN source activate hecaton_py2 && \
	git clone --recursive https://github.com/hall-lab/speedseq && \
	cd speedseq && \
	make align && \
	make sv && \
	make config && \
	cd .. && \
	source deactivate

# Setup hecaton pipeline
# RUN git clone https://git.wur.nl/wijfj001/hecaton.git && cd hecaton && \
#	git checkout e85bba0c && cd .. && \
RUN git clone https://git.wur.nl/bioinformatics/hecaton.git && \
	echo "f82885e8" && \
	chmod +x hecaton/scripts/collapse/* && \
	chmod +x hecaton/scripts/convert/* && \
	chmod +x hecaton/scripts/filter/* && \
	chmod +x hecaton/scripts/genotype/* && \ 
	chmod +x hecaton/scripts/gridss/* && \
	chmod +x hecaton/scripts/intersect/* && \
	chmod +x hecaton/scripts/predict/* && \
	chmod +x hecaton/scripts/process/*  

# Set paths
ENV PATH /opt/conda/envs/hecaton/bin:/speedseq/bin:/hecaton/scripts/collapse:/hecaton/scripts/convert:/hecaton/scripts/filter:/hecaton/scripts/genotype:/hecaton/scripts/gridss:/hecaton/scripts/intersect:/hecaton/scripts/predict:/hecaton/scripts/process:$PATH
ENV PYTHONPATH /hecaton/scripts

# Clean up
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

