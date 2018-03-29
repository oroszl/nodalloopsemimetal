FROM jupyter/scipy-notebook:cf6258237ff9
USER root
RUN apt-get update
RUN apt-get install -y cmake
USER jovyan
RUN conda install numpy scipy matplotlib
RUN conda install ipyvolume 
RUN conda install tqdm
RUN pip install pybinding --user --upgrade
