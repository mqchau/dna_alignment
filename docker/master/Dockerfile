FROM ubuntu:14.04

RUN apt-get install python python-dev pip libpq-dev git wget
RUN pip install psycopg2 sqlalchemy numpy boto3

# copy the aws credentials over
RUN mkdir /root/.aws
ADD dont_commit_aws_credentials /root/.aws/credentials
ADD dont_commit_aws_config /root/.aws/config

# clone the git repo
RUN git clone https://github.com/mqchau/dna_alignment /tmp

WORKDIR /tmp

# download the ref and reads file
RUN mkdir -p /tmp/dataset/ugrad
RUN wget -O /tmp/dataset/ugrad/ref_ugrad.txt https://s3-us-west-2.amazonaws.com/cm124-hw/ref_ugrad.txt
RUN wget -O /tmp/dataset/ugrad/reads_ugrad.txt https://s3-us-west-2.amazonaws.com/cm124-hw/reads_ugrad.txt

# set up env for RDS instance
ENV RDS_HOST cm124-hw-db.cujwspoan2k0.us-west-2.rds.amazonaws.com
ENV RDS_PORT 5432
ENV RDS_DATABASE bioinfo

# run worker loop
CMD ["python", "master_loop.py"]


