# Dockerfile for Fragment Explorer
# author: Micha Birklbauer
# version: 1.0.0

FROM ubuntu:22.04

LABEL maintainer="micha.birklbauer@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    git \
    python3-distutils \
    python3-lxml \
    python3-pip

RUN git clone https://github.com/veitveit/eubic23_internal_ions.git
WORKDIR eubic23_internal_ions

RUN pip3 install -r requirements.txt
RUN pip3 install -r fragannot_requirements.txt

CMD  ["streamlit", "run", "streamlit_app.py"]
