FROM amazonlinux:2

ARG AWS_REGION

RUN yum -y update
RUN yum -y install python37 tar git unzip sudo jq
RUN curl -O https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py --user

RUN curl -OLS https://github.com/BenLangmead/bowtie2/releases/download/v2.4.2/bowtie2-2.4.2-linux-x86_64.zip
RUN unzip bowtie2-2.4.2-linux-x86_64.zip

COPY bowtie2-2.4.2-linux-x86_64/* /usr/local/bin/

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN sudo ./aws/install

RUN aws configure set default.region us-east-1
RUN git clone https://github.com/sternberglab/Illumina-pipeline
RUN pip3 install -r Illumina-pipeline/requirements.txt

CMD cd ./Illumina-pipeline && /bin/bash ./ingestion.sh