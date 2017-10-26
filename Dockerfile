from centos:6

# Required packages
RUN yum install -y wget bzip2 libXext libSM libXrender libpng libgomp

RUN mkdir -p /src /app /NGSData

# Install
ADD . /src
WORKDIR /src
RUN sed 's%/path/to%%' ngs_mapper/config.yaml.default > ngs_mapper/config.yaml
RUN ./install.sh /app/miniconda

# Set matplotlib backend to Agg so it doesn't require DISPLAY
WORKDIR /app
RUN ls -l /app
RUN sed -i 's/Qt4Agg/Agg/' miniconda/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc

ENV PATH /app/miniconda/bin:$PATH

# Cleanup to make image smaller
RUN yum clean all
RUN rm -rf /src
