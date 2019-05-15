FROM debian:latest
MAINTAINER Yifan Zhang

# Update OS
RUN apt-get update -y

# Install Python other things
#RUN apt-get install -y python-pip python-dev build-essential
RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    vim \
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

RUN mkdir -p log
RUN touch log/bartweb_backend.log
RUN chown -R www-data:www-data log
RUN chmod -R 775 log

#ENTRYPOINT ["sh"]
#CMD ["run.sh"]