FROM jupyter/r-notebook

LABEL maintainer="illumination-k <illumination.k.27@gamil.com>"

USER root

LABEL illumination-k <illumination.k.27@gmail.com>

RUN apt-get update --fix-missing && \
    apt-get install -y locales && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8

ENV LANG=en_US.UTF-8
ENV LC_CTYPE=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

RUN conda install \
    bioconductor-edger \
    bioconductor-deseq2 \
    bioconductor-goseq \
    bioconductor-topgo \
    bioconductor-rgraphviz \
    bioconductor-genefilter \
    bioconductor-annotationdbi \
    bioconductor-clusterprofiler \
    bioconductor-sva \
    r-tidyverse \
    r-wgcna \
    pygraphviz \
    -c bioconda -c r -c conda-forge

RUN mkdir -p /workspace
WORKDIR /workspace

ADD ./requirements.txt /workspace

RUN pip install -U pip && \
    pip install -r requirements.txt

RUN jupyter labextension install @lckr/jupyterlab_variableinspector \
                                 @jupyterlab/git \
                                 @jupyterlab/toc \
                                 @jupyter-widgets/jupyterlab-manager \
                                 jupyterlab-plotly \
                                 plotlywidget \
                                 --no-build && \
    pip install jupyterlab-git && \
    jupyter serverextension enable --py jupyterlab_git && \
    jupyter lab build

RUN pip install jupyterthemes

# black background
RUN mkdir -p /home/jovyan/.jupyter/lab/user-settings/@jupyterlab/apputils-extension
RUN echo '{"theme":"JupyterLab Dark"}' > \
  /home/jovyan/.jupyter/lab/user-settings/@jupyterlab/apputils-extension/themes.jupyterlab-settings

EXPOSE 8888