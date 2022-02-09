FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

ENV HMMER_VERSION='3.3.2'
ENV dbCAN_VERSION='10'
#ENV ETE3_VERSION='3.0.0b35'
ENV ETE3_VERSION='3.1.1'
#ENV ETE3_VERSION='3.1.2'

# Update packages
RUN apt-get update

# Install ETE3 dependencies
#RUN apt-get update && \
#    apt-get -y install xvfb python-qt4 && \
#    pip install ete3==${ETE3_VERSION}
RUN apt-get -y install xvfb python-qt4 python-numpy python-lxml python-six
    
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

# install ETE3
#
WORKDIR /kb/module
RUN \
  git clone -b ${ETE3_VERSION} https://github.com/etetoolkit/ete


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
  curl https://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V${dbCAN_VERSION}.txt > data/dbCAN/dbCAN-v${dbCAN_VERSION}/dbCAN-fam-HMMs-v${dbCAN_VERSION}.txt && \
  cd ete && \
  python setup.py install


# Start up
#
WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
