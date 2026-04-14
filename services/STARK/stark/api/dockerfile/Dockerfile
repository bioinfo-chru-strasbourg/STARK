FROM almalinux:9

# Update
#RUN yum update -y && yum install -y epel-release && yum install -y python3 python-pip3 docker git make gcc g++ && yum clean all
RUN yum update -y && yum install -y epel-release && yum install -y python3-pip wget git make gcc g++ && yum clean all

# Install app dependencies
COPY requirements.txt /src/requirements.txt
RUN pip3 install -r /src/requirements.txt


ENV TOOL_NAME="docker"
ENV TOOL_VERSION="29.2.1"
ENV TOOL_TARBALL="docker-$TOOL_VERSION.tgz"
ENV TOOL_SOURCE_EXTERNAL="https://download.docker.com/linux/static/stable/x86_64/$TOOL_TARBALL"
ENV PATH=$PATH:/tools/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS
ENV DOCKER_HOST=unix:///var/run/docker.sock

# TOOL INSTALLATION
RUN echo "#[INFO] SYSTEM Docker installation '$TOOL_NAME:$TOOL_VERSION'"

RUN mkdir -p /tools/$TOOL_NAME/$TOOL_VERSION/bin && \
    echo wget -q $TOOL_SOURCE_EXTERNAL -O $TOOL_TARBALL && \
    wget -q $TOOL_SOURCE_EXTERNAL -O $TOOL_TARBALL && \
    tar -xvf $TOOL_TARBALL -C /tools/$TOOL_NAME/$TOOL_VERSION/bin --strip-components=1 && \
    ln -s /tools/$TOOL_NAME/$TOOL_VERSION/bin/docker /usr/local/bin/docker && \
    rm -rf $TOOL_TARBALL;


# Task Spooler
RUN git clone https://github.com/thomaspreece/task-spooler.git && \
    cd task-spooler && \
    make && \
    make install && \
    mkdir /ts-tmp && \
    rm -rf task-spooler ;


# Bundle app source
COPY . /app

WORKDIR /app


# YUM clean
RUN yum erase -y git make gcc g++ && yum clean all ;


EXPOSE  8000
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
