FROM python:3.6
ENV PYTHONUNBUFFERED 1
RUN mkdir /clust_app
WORKDIR /clust_app
ADD . /clust_app/
RUN pip install -r requirements.txt
