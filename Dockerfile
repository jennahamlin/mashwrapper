FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
	apt-get install -y wget python3.6 python3-distutils python3-pip python3-pandas python3-tabulate && \
	apt-get clean


RUN wget https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar && \
    tar -xvf mash-Linux64-v2.0.tar && \
    rm -rf mash-Linux64-v2.0.tar

ENV PATH="${PATH}:/mash-Linux64-v2.0" \
    LC_ALL=C
