from centos:6

# Required packages
RUN yum install -y wget bzip2 libXext libSM libXrender libpng libgomp

RUN mkdir -p /app /NGSData
ADD . /app
WORKDIR /app
RUN ./install.sh miniconda

# Set matplotlib backend to Agg so it doesn't require DISPLAY
RUN sed -i 's/Qt4Agg/Agg/' miniconda/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc

ENV PATH /app/miniconda/bin:$PATH
