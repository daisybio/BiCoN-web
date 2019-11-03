# pull official python base image
FROM python:3.6
# Unbuffered stdin / out for faster duping of logs
ENV PYTHONUNBUFFERED 1

# create the bigants user with home directory
RUN addgroup --system bigants && adduser --home /home/bigants --system bigants --ingroup bigants

# create the appropriate directories
ENV HOME=/home/bigants
ENV APP_HOME=/home/bigants/web
RUN mkdir $APP_HOME


# set work directory
WORKDIR $APP_HOME

# install dependencies
COPY requirements/requirements.txt $APP_HOME/requirements/
RUN pip install -r requirements/requirements.txt

## TEMPORARILY START
#COPY _dev/BiGAnts-PyPI-package $APP_HOME/_dev/BiGAnts-PyPI-package
#WORKDIR $APP_HOME/_dev/BiGAnts-PyPI-package
#RUN python setup.py install
#WORKDIR $APP_HOME
## TEMPORARILY END

# copy project
COPY . $APP_HOME

# chown all the files to the bigants user
RUN chown -R bigants:bigants $APP_HOME

# change to the bigants user
USER bigants

# run entrypoint.prod.sh
#ENTRYPOINT ["/home/app/web/entrypoint.prod.sh"]