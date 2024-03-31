FROM python:3.10.14-slim-bullseye
RUN apt-get update \
    && apt-get install -y curl vim less procps git openssh-server sudo && \
    adduser --disabled-password --gecos '' python && \
    gpasswd -a python sudo && passwd -d python && \
    curl https://packages.microsoft.com/keys/microsoft.asc | tee /etc/apt/trusted.gpg.d/microsoft.asc && \
    curl https://packages.microsoft.com/config/debian/11/prod.list | tee /etc/apt/sources.list.d/mssql-release.list

USER python
WORKDIR /workspaces
#RUN pip install --upgrade pip && pip install --user --no-cache-dir -r requirements.txt
CMD /bin/sh -c "while sleep 1000; do :; done"