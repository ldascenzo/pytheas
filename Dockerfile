# Start from Ubuntu 20.04
FROM ubuntu:20.04

MAINTAINER Luigi <dascenzoluigi@gmail.com>

RUN apt-get update && apt-get install -y \
    software-properties-common
RUN add-apt-repository universe
RUN apt-get update && apt-get install -y \
	wget \
	python3 \
	python3-pip
RUN pip3 install Biopython pandas lxml numpy matplotlib
RUN pip3 install xlrd==1.2.0 pyteomics==4.1.2
RUN wget  -P ~ "https://extras.wxpython.org/wxPython4/extras/linux/gtk3/ubuntu-20.04/wxPython-4.1.0-cp38-cp38-linux_x86_64.whl"
RUN pip3 install ~/wxPython-4.1.0-cp38-cp38-linux_x86_64.whl
RUN pip3 install gooey

#RUN mkdir -p Pytheas

VOLUME ./Digest