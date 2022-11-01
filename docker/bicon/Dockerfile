# pull official python base image
FROM python:3.10
# Unbuffered stdin / out for faster duping of logs
ENV PYTHONUNBUFFERED 1

# create the bicon user with home directory
RUN addgroup --system bicon && adduser --home /home/bicon --system bicon --ingroup bicon

# create the appropriate directories
ENV HOME=/home/bicon
ENV APP_HOME=/home/bicon/web
RUN mkdir $APP_HOME


# set work directory
WORKDIR $APP_HOME

# install dependencies
COPY requirements/requirements.txt $APP_HOME/requirements/
RUN pip install -r requirements/requirements.txt

## TEMPORARILY START -- Install BiGAnts package from local source
#COPY _dev/BiGAnts-PyPI-package $APP_HOME/_dev/BiGAnts-PyPI-package
#WORKDIR $APP_HOME/_dev/BiGAnts-PyPI-package
#RUN python setup.py install
#WORKDIR $APP_HOME
## TEMPORARILY END

# copy project
COPY . $APP_HOME

# chown all the files to the bicon user
RUN chown -R bicon:bicon $APP_HOME

# change to the bicon user
USER bicon

# run entrypoint.prod.sh
#ENTRYPOINT ["/home/app/web/entrypoint.prod.sh"]