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
RUN mkdir -p /BART_backend
COPY ./ /BART_backend/
RUN pip3 install -r /BART_backend/requirements.txt
WORKDIR /BART_backend



#ENTRYPOINT ["sh"]
#CMD ["run.sh"]