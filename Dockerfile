FROM python:3.6
ENV PYTHONUNBUFFERED 1
RUN mkdir /testproject_2
WORKDIR /testproject_2
ADD . /testproject_2/
RUN pip install -r requirements.txt
