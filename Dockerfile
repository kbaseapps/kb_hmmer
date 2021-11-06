FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

ENV HMMER_VERSION='3.3.2'
ENV dbCAN_VERSION='10'

# Update packages
RUN apt-get update

# add packages
#RUN apt install -y build-essential

# Here we install a python coverage tool and an
# https library that is out of date in the base image.
RUN pip install coverage


# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

# Install HMMER
#
WORKDIR /kb/module
RUN rm -f /usr/bin/hmm*
RUN \
  curl http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz > hmmer-${HMMER_VERSION}.tar.gz && \
  tar xfz hmmer-${HMMER_VERSION}.tar.gz && \
  ln -s hmmer-${HMMER_VERSION} hmmer && \
  rm -f hmmer-${HMMER_VERSION}.tar.gz
WORKDIR /kb/module/hmmer
RUN \
  ./configure --prefix /kb/module/hmmer && \
  make && \
  make install


# Install dbCAN HMM data
#
WORKDIR /kb/module
RUN \
  curl https://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V${dbCAN_VERSION}.txt > data/dbCAN/dbCAN-v${dbCAN_VERSION}/dbCAN-fam-HMMs-v${dbCAN_VERSION}.txt


# Start up
#
WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
