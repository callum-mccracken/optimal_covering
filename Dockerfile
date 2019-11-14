FROM continuumio/miniconda3

LABEL Name=optimal_coverings Version=0.0.1
EXPOSE 3000

ADD environment.yml /tmp/environment.yml

RUN conda env create -f environment.yml
CMD /bin/bash -c "source activate optcov"
