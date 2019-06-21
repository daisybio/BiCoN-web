FROM python:3.6
ENV PYTHONUNBUFFERED 1
RUN mkdir /clust_app
WORKDIR /clust_app
ADD . /clust_app/
RUN pip install -r requirements.txt
ADD crontab /etc/cron.d/hello-cron
RUN chmod 0644 /etc/cron.d/hello-cron
RUN touch /var/log/cron.log
RUN apt-get update
RUN apt-get -y install cron
CMD cron && tail -f /var/log/cron.log
RUN cp crontab /etc/crontab
