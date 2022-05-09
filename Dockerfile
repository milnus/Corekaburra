FROM python:3.9.7-buster
WORKDIR /Corekaburra
COPY . .

# Install the python package (and executable)
RUN pip3 install .

# Override some of the dependencies with the hard-coded versions
RUN pip3 install -r requirements-dev.txt
