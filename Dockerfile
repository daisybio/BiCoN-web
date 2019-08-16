FROM python:3.6
ENV PYTHONUNBUFFERED 1
RUN mkdir /clust_app
WORKDIR /clust_app
COPY . /clust_app/
RUN pip install -r requirements.txt
COPY crontab /etc/cron.d/hello-cron
RUN chmod 0644 /etc/cron.d/hello-cron
RUN touch /var/log/cron.log
RUN apt-get update
RUN apt-get -y install cron
RUN cp crontab /etc/crontab
COPY crontab /etc/cron.d/cool-task
RUN chmod 0644 /etc/cron.d/cool-task
CMD cron && tail -f /var/log/cron.log && crontab crontab && service cron start
