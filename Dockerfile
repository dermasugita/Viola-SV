FROM python:3.10.14-slim-bullseye
RUN apt-get update \
    && apt-get install -y curl vim less procps git openssh-server sudo build-essential && \
    adduser --disabled-password --gecos '' python && \
    gpasswd -a python sudo && passwd -d python && \
    curl https://packages.microsoft.com/keys/microsoft.asc | tee /etc/apt/trusted.gpg.d/microsoft.asc && \
    curl https://packages.microsoft.com/config/debian/11/prod.list | tee /etc/apt/sources.list.d/mssql-release.list

USER python
WORKDIR /workspaces
COPY setup.py .
COPY pyproject.toml .
COPY README.rst .
COPY ./src/viola/_version.py ./src/viola/_version.py
RUN sudo chown -R python:python ./src && pip install --upgrade pip setuptools && pip install --editable .
#RUN pip install --upgrade pip && pip install --user --no-cache-dir -r requirements.txt
CMD /bin/sh -c "while sleep 1000; do :; done"