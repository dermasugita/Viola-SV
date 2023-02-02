FROM mcr.microsoft.com/devcontainers/python:0-3.9
USER vscode
ENV HOME /home/vscode
ENV PYENV_ROOT $HOME/.pyenv
ENV PATH $PYENV_ROOT/bin:$PATH
RUN git clone https://github.com/pyenv/pyenv.git ~/.pyenv && \
    echo 'eval "$(pyenv init -)"' >> ~/.bashrc
RUN pyenv install 3.6.14 3.7.15 3.8.16 3.9.16 3.10.9 3.11.1
RUN pyenv global 3.6.14 3.7.15 3.8.16 3.9.16 3.10.9 3.11.1
RUN pip install tox