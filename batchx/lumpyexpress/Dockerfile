# FROM gcc:latest
# WORKDIR /batchx
# # samtools
# RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2
# RUN tar -xjvf samtools-1.15.tar.bz2 && samtools-1.15/configure && make install -C samtools-1.15
# # samblaster
# RUN git clone https://github.com/GregoryFaust/samblaster.git
# RUN make -C samblaster && cp samblaster/samblaster /usr/local/bin/.
# # lumpy
# RUN git clone --recursive https://github.com/arq5x/lumpy-sv.git
# RUN make -C lumpy-sv && cp lumpy-sv/bin/* /usr/local/bin/.
# # other dependencies
# RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python get-pip.py
# RUN pip install pysam numpy jq
# RUN apt-get update && apt-get install -y gawk tabix && apt -y autoremove
# RUN apt install -y bsdmainutils
# RUN rm /opt/conda/bin/lumpyexpress
# COPY lumpyexpress /opt/conda/bin/lumpyexpress
# COPY lumpyexpress.sh /batchx/lumpyexpress.sh
# RUN chmod -R 777 /batchx
# ENTRYPOINT /batchx/lumpyexpress.py
# LABEL io.batchx.manifest=08
# COPY manifest /batchx/manifest/

FROM continuumio/miniconda2
WORKDIR /batchx
# samblaster
RUN conda install -c bioconda samblaster
RUN conda install -c bioconda lumpy-sv
RUN conda install -c conda-forge jq
RUN apt install -y bsdmainutils
RUN rm /opt/conda/bin/lumpyexpress
COPY lumpyexpress /opt/conda/bin/lumpyexpress
COPY lumpyexpress.sh /batchx/lumpyexpress.sh
RUN chmod -R 777 /batchx
ENTRYPOINT bash /batchx/lumpyexpress.sh
LABEL io.batchx.manifest=10
COPY manifest /batchx/manifest/

# FROM --platform=linux/amd64 gcc:latest
# WORKDIR /batchx
# RUN apt update && apt-get -y install python2 build-essential wget git libcurl4-gnutls-dev libxml2-dev libssl-dev autoconf && apt-get autoremove -y
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda
# ENV PATH=$CONDA_DIR/bin:$PATH
# RUN conda install -c conda-forge jq
# # samtools
# # RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2
# # RUN tar -xjvf samtools-1.15.tar.bz2 && samtools-1.15/configure && make install -C samtools-1.15
# # samblaster
# RUN conda install -c bioconda samblaster
# RUN apt install -y bsdmainutils
# # lumpy
# RUN conda install -c bioconda lumpy-sv
# # sambamba
# # RUN conda install -c bioconda sambamba
# # other dependencies
# # RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python get-pip.py
# # RUN pip install pysam numpy
# # RUN apt-get update && apt-get install -y gawk tabix && apt -y autoremove
# # lumpy
# # RUN conda install -c bioconda lumpy-sv
# COPY lumpyexpress /opt/conda/bin/lumpyexpress
# COPY lumpyexpress.sh /batchx/lumpyexpress.sh
# RUN chmod -R 777 /batchx
# ENTRYPOINT sh /batchx/lumpyexpress.sh
# LABEL io.batchx.manifest=08
# COPY manifest /batchx/manifest/

# FROM jbwebster/lumpy_docker
# WORKDIR /batchx
# RUN apt update && apt-get -y install build-essential wget libcurl4-gnutls-dev bsdmainutils && apt-get autoremove -y
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#      /bin/bash ~/miniconda.sh -b -p /opt/conda
# ENV PATH=$CONDA_DIR/bin:$PATH
# RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py && python get-pip.py
# RUN pip install pysam numpy pysmoove hexdump
# RUN conda install -c conda-forge jq
# RUN apt install -y bsdmainutils
# COPY lumpyexpress /opt/conda/bin/lumpyexpress
# COPY lumpyexpress.sh /batchx/lumpyexpress.sh
# RUN conda create --name py2 python=2.7
# RUN chmod -R 777 /batchx
# ENTRYPOINT sh /batchx/lumpyexpress.sh
# LABEL io.batchx.manifest=08
# COPY manifest /batchx/manifest/

# libcurl4-gnutls-dev libxml2-dev libssl-dev

# FROM --platform=linux/amd64 python:2.7 
# WORKDIR /batchx
# # samblaster
# RUN apt update && apt-get -y install build-essential wget && apt-get autoremove -y
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda
# ENV PATH=$CONDA_DIR/bin:$PATH
# RUN conda install -c conda-forge jq
# RUN apt install -y bsdmainutils
# RUN conda install -c bioconda samblaster
# RUN apt install python2
# RUN conda install -c bioconda lumpy-sv
# RUN rm /opt/conda/bin/lumpyexpress
# COPY lumpyexpress /opt/conda/bin/lumpyexpress
# COPY lumpyexpress.sh /batchx/lumpyexpress.sh
# RUN chmod -R 777 /batchx
# ENTRYPOINT sh /batchx/lumpyexpress.sh
# LABEL io.batchx.manifest=08
# COPY manifest /batchx/manifest/

# git libcurl4-gnutls-dev libxml2-dev libssl-dev autoconf 

