# pull official python base image
FROM python:3.6
# Unbuffered stdin / out for faster duping of logs
ENV PYTHONUNBUFFERED 1

# create directory for the app user
RUN mkdir -p /home/bigants

# create the bigants user
RUN addgroup -S bigants && adduser -S bigants -G bigants

# create the appropriate directories
ENV HOME=/home/bigants
ENV APP_HOME=/home/bigants/web
RUN mkdir $APP_HOME


# set work directory
WORKDIR $APP_HOME

# install dependencies
COPY requirements/requirements.txt /django/requirements/
RUN pip install -r requirements/requirements.txt

# copy entrypoint-prod.sh
COPY ./entrypoint.prod.sh $APP_HOME

# copy project
COPY . $APP_HOME

# chown all the files to the bigants user
RUN chown -R bigants:bigants $APP_HOME

# change to the bigants user
USER bigants

# run entrypoint.prod.sh
#ENTRYPOINT ["/home/app/web/entrypoint.prod.sh"]