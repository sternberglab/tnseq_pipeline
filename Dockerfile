FROM amazonlinux:1

ARG AWS_REGION

COPY . /root/pipeline

RUN yum -t update
RUN yum -y install python36 tar git
RUN rm /usr/bin/python
RUN ls -s /etc/alternatives/python3 /usr/bin/python

RUN aws configure set default.region $AWS_REGION
RUN git clone https://github.com/sternberglab/Illumina-pipeline -y
RUN pip install -r Illumina-pipeline/requirements.txt

CMD cd /root/pipeline && /bin/bash ./ingestion.sh