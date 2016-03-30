FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.



# Update Transform (should go away eventually)
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf transform && \ 
  git clone https://github.com/kbase/transform -b develop

# setup the transform, but ignore errors because sample data cannot be found!
RUN \
  . /kb/dev_container/user-env.sh; \
  cd /kb/dev_container/modules/transform/t/demo; \
  python setup.py; \
  exit 0;



# RUN apt-get update

# -----------------------------------------

# Install SDK Module
#
RUN mkdir -p /kb/module
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
WORKDIR /kb/module
RUN make


# Install HMMER
#
WORKDIR /kb/module
RUN \
  curl http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz > hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
  tar xfz hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
  ln -s hmmer-3.1b2-linux-intel-x86_64 hmmer && \
  rm -f hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
  cd hmmer && \
  ./configure && \
  make


WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
