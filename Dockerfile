FROM jupyter/scipy-notebook:cf6258237ff9
RUN conda install numpy scipy matplotlib
RUN conda install ipyvolume 
RUN conda install tqdm
RUN pip install pybinding --user --upgrade