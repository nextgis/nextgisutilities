FROM ubuntu:bionic

ARG DEBIAN_FRONTEND=noninteractive

ENV LOGIN=user
ENV PASSWORD=secret1

RUN echo "machine rm.nextgis.com\nlogin $LOGIN\npassword $PASSWORD\n" > /etc/apt/auth.conf.d/rm.conf
RUN apt-get update -y && \
    apt-get -y install --no-install-recommends --yes \
    apt-transport-https ca-certificates curl gnupg && \
    echo "deb https://rm.nextgis.com/api/repo/11/deb bionic main" | tee -a /etc/apt/sources.list && \
    curl -s -L https://rm.nextgis.com/api/repo/11/deb/key.gpg | apt-key add - && \
    echo "deb https://rm.nextgis.com/api/repo/12/deb bionic main" | tee -a /etc/apt/sources.list && \
    curl -s -L -u "$LOGIN:$PASSWORD" https://rm.nextgis.com/api/repo/12/deb/key.gpg | apt-key add - && \
    apt-get update -y && \
    apt-get -y install --no-install-recommends --yes bash nextgisutilities-bin gdal-bin

RUN mkdir /data
WORKDIR /data
