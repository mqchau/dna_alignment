FROM ubuntu:14.04

RUN apt-get install python python-dev pip libpq-dev git
RUN pip install psycopg2 sqlalchemy numpy boto3

# copy the aws credentials over
RUN mkdir /root/.aws
ADD dont_commit_aws_credentials /root/.aws/credentials
ADD dont_commit_aws_config /root/.aws/config

# clone the git repo
RUN git clone https://github.com/mqchau/dna_alignment /tmp

WORKDIR /tmp

# set up env for RDS instance
ENV RDS_HOST cm124-hw.cujwspoan2k0.us-west-2.rds.amazonaws.com
ENV RDS_PORT 5432
ENV RDS_DATABASE bioinfo

# run worker loop
CMD ["python", "worker_loop.py"]


