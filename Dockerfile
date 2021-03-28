FROM amazonlinux:2

ARG AWS_REGION

RUN yum -y update
RUN yum -y install python37 tar git unzip sudo
RUN curl -O https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py --user

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN sudo ./aws/install

RUN aws configure set default.region us-east-1
RUN git clone https://github.com/sternberglab/Illumina-pipeline
RUN pip3 install -r Illumina-pipeline/requirements.txt

CMD cd /root/Illumina-pipeline && /bin/bash ./ingestion.sh