FROM almalinux:9

# Update
#RUN yum update -y && yum install -y epel-release && yum install -y python3 python-pip3 docker git make gcc g++ && yum clean all
RUN yum update -y && yum install -y epel-release && yum install -y python3-pip wget git make gcc g++ && yum clean all

# Install app dependencies
COPY requirements.txt /src/requirements.txt
RUN pip3 install -r /src/requirements.txt

# Docker installation
ENV TOOL_NAME="docker"
ENV TOOL_VERSION="29.2.1"
ENV TOOL_TARBALL="docker-$TOOL_VERSION.tgz"
ENV TOOL_SOURCE_EXTERNAL="https://download.docker.com/linux/static/stable/x86_64/$TOOL_TARBALL"
ENV PATH=$PATH:/tools/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS
ENV DOCKER_HOST=unix:///var/run/docker.sock

# TOOL INSTALLATION
RUN echo "#[INFO] SYSTEM Docker installation '$TOOL_NAME:$TOOL_VERSION'" && \
    mkdir -p /tools/$TOOL_NAME/$TOOL_VERSION/bin && \
    wget -q $TOOL_SOURCE_EXTERNAL -O $TOOL_TARBALL && \
    tar -xvf $TOOL_TARBALL -C /tools/$TOOL_NAME/$TOOL_VERSION/bin --strip-components=1 && \
    ln -s /tools/$TOOL_NAME/$TOOL_VERSION/bin/docker /usr/local/bin/docker && \
    rm -rf $TOOL_TARBALL;


# Docker compose installation
ENV TOOL_NAME="docker-compose"
ENV TOOL_VERSION="v5.1.3"
ENV TOOL_TARBALL="docker-compose-$TOOL_VERSION.tgz"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/docker/compose/releases/download/$TOOL_VERSION"
ENV PATH=$PATH:/tools/$TOOL_NAME/$TOOL_VERSION/bin
# TOOL PARAMETERS

# TOOL INSTALLATION
RUN echo "#[INFO] SYSTEM Docker installation '$TOOL_NAME:$TOOL_VERSION'" && \
    mkdir -p /tools/$TOOL_NAME/$TOOL_VERSION/bin && \
    wget -q $TOOL_SOURCE_EXTERNAL/docker-compose-$(uname -s)-$(uname -m) -O /tools/$TOOL_NAME/$TOOL_VERSION/bin/docker-compose && \
    chmod +x /tools/$TOOL_NAME/$TOOL_VERSION/bin/docker-compose && \
    ln -s /tools/$TOOL_NAME/$TOOL_VERSION/bin/docker-compose /usr/local/bin/docker-compose;


# Task Spooler
ENV TOOL_NAME="task-spooler"
ENV TOOL_VERSION="v0.5.4"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/thomaspreece/task-spooler.git"

# TOOL INSTALLATION
RUN echo "#[INFO] SYSTEM Task Spooler installation '$TOOL_NAME:$TOOL_VERSION'" && \
    #git clone $TOOL_SOURCE_EXTERNAL && \
    git clone --branch $TOOL_VERSION --depth 1 $TOOL_SOURCE_EXTERNAL && \
    cd task-spooler && \
    make && \
    make install && \
    mkdir /ts-tmp && \
    rm -rf task-spooler ;


# Bundle app source
COPY app /app

WORKDIR /app


# YUM clean
RUN yum erase -y git make gcc g++ && yum clean all ;


EXPOSE  8000
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
