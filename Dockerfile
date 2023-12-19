FROM python:3.8.12

MAINTAINER Firas Hamood "firas.hamood@tum.de"

LABEL website=https://gitlab.lrz.de/proteomics/clustering_transfer_tool
LABEL description="SIMSI-Transfer"
LABEL tags="mass spectrometry tmt isobaric labeling"
LABEL documentation=https://gitlab.lrz.de/proteomics/clustering_transfer_tool

# Tell docker that we don't want to be bothered with questions
ARG DEBIAN_FRONTEND=noninteractive

# for mono installation
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/debian stable-buster main" | tee /etc/apt/sources.list.d/mono-official-stable.list

RUN apt-get update && apt-get install -y \
        mono-devel \
        zip \
    && rm -rf /var/lib/apt/lists/*

# set root directory
ENV HOME /root
WORKDIR /root

RUN pip install poetry==1.5.1
# poetry uses virtualenvs by default -> we want global installation
RUN poetry config virtualenvs.create false
ADD pyproject.toml /root/pyproject.toml
ADD poetry.lock /root/poetry.lock
RUN poetry install


# install maracluster
RUN ZIP=ubuntu.tar.gz && \
    wget https://github.com/statisticalbiotechnology/maracluster/releases/download/rel-1-01-1/$ZIP -O /tmp/$ZIP && \
    tar xvzf /tmp/$ZIP && \
    chmod -R 755 /tmp/* && \
    dpkg -i maracluster-v1-01-linux-amd64.deb && \
    rm /tmp/$ZIP

RUN pip install memory_profiler

ADD simsi_transfer/ /root/simsi_transfer
