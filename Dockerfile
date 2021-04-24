FROM continuumio/miniconda

RUN conda install python=3.8
RUN apt-get update -y
RUN apt-get -y install gcc make libuv1-dev

##Grab BLAST
WORKDIR /source
COPY {pkg} /source/{pkg}/
COPY setup.py /source/.
RUN pip3 install .
