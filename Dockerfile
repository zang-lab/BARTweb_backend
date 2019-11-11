FROM debian:latest
MAINTAINER Yifan Zhang

# Update OS
RUN apt-get update -y

# Install Python other things
RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    vim \
    jq \
 && apt-get clean \
 && apt-get autoremove \
 && rm -rf /var/lib/apt/lists/*

# && do not continue if any fails

# ADD . /app
RUN mkdir -p /BARTweb
COPY ./ /BARTweb/
RUN pip3 install -r /BARTweb/requirements.txt
WORKDIR /BARTweb

# give docker permission to write log files
RUN mkdir -p usercase/log
RUN touch usercase/log/bartweb_backend.log
RUN chown -R www-data:www-data usercase/log
RUN chmod -R 775 usercase/log

ENV HDF5_USE_FILE_LOCKING="FALSE"

#ENTRYPOINT ["sh"]
#CMD ["run.sh"]
#CMD ["check.sh"]
