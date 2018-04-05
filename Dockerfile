FROM jupyter/scipy-notebook:cf6258237ff9
USER root
RUN apt-get update
RUN apt-get install -y cmake
USER jovyan
RUN conda install numpy scipy 
RUN conda install matplotlib=2.1.1
RUN conda install ipyvolume 
RUN conda install tqdm
RUN conda 
RUN pip install pybinding --user --upgrade
RUN git clone https://github.com/oroszl/nodalloopsemimetal
