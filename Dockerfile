FROM ubuntu:18.04

# Upgrade installed packages
RUN apt-get update && apt-get upgrade -y && apt-get clean

# Python package management and basic dependencies
RUN apt-get install -y curl wget python3.7 python3.7-dev python3.7-distutils

# Register the version in alternatives
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1

# Set python 3 as the default python
RUN update-alternatives --set python /usr/bin/python3.7

# Upgrade pip to latest version
RUN curl -s https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py --force-reinstall && \
    rm get-pip.py

RUN pip3 install pandas tabulate

RUN wget https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar && \
    tar -xvf mash-Linux64-v2.0.tar && \
    rm -rf mash-Linux64-v2.0.tar

ENV PATH="${PATH}:/mash-Linux64-v2.0"

##From: https://stackoverflow.com/questions/56135497/can-i-install-python-3-7-in-ubuntu-18-04-without-having-python-3-6-in-the-system
