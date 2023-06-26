FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade && apt-get install -y build-essential wget curl unzip cmake python3-pip mono-complete python-is-python3

COPY . /root/protgraph_identification
WORKDIR /root/protgraph_identification

RUN ./compile_and_setup_depencies.sh
ENV PATH="${PATH}://root/protgraph_identification/bin"
