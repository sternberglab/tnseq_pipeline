FROM amazonlinux:2

ARG AWS_REGION

COPY . /root/pipeline

RUN yum -y update
RUN yum -y install python36 tar git
RUN rm /usr/bin/python
RUN ln -s /etc/alternatives/python3 /usr/bin/python

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN sudo ./aws/install

RUN aws configure set default.region $AWS_REGION
RUN git clone https://github.com/sternberglab/Illumina-pipeline -y
RUN pip install -r Illumina-pipeline/requirements.txt

CMD cd /root/pipeline && /bin/bash ./ingestion.sh