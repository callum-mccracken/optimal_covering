FROM continuumio/miniconda3

LABEL Name=optimal-coverings Version=0.0.1
EXPOSE 3000

ADD environment.yml environment.yml

RUN conda env create -f environment.yml
CMD /bin/bash -c "conda activate optcov"
