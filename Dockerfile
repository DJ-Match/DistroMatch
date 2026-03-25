FROM ubuntu:22.04
RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker
RUN echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker
RUN DEBIAN_FRONTEND=noninteractive \
  apt-get update \
  && apt-get install -y python3 \
  && apt-get install -y pip \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app
RUN useradd -ms /bin/bash apprunner
RUN chown apprunner:apprunner /app
USER apprunner


RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir numpy pandas matplotlib seaborn scipy PyYAML networkit

COPY . /app
