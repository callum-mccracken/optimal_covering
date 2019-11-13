FROM continuumio/miniconda3

LABEL Name=optimal_coverings Version=0.0.1
EXPOSE 3000

WORKDIR /app
ADD . /app

# Using miniconda (make sure to replace 'myenv' w/ your environment name):
RUN conda env create -f environment.yml
CMD /bin/bash -c "source activate myenv && python3 -m optimal_coverings"
