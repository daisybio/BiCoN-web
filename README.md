# 1. Getting Started
BiGAnts-web comes in an easy to use `docker-compose` format.

To deploy BiGAnts-web in a production environment, please ensure you have a working docker (and docker-compose) environment running. 
You find a manual for install docker on your operating system here.

Once docker is installed and running, execute the following commands:

## Running BiGAnts-web for the first time
To install BiGAnts-web (running it for the first time) follow these steps:
```shell script
# First clone this repository and change into the created directory
git clone https://github.com/biomedbigdata/BiGAnts-web.git && cd BiGAnts-web

# Change the branch to `new_structure`
git checkout new_structure

# Change the password and username for the databases (maybe leave this out, Markus?)

# Deploy and build the containers
docker-compose up -d --build

# Apply migrations to the database
docker-compose exec web python manage.py migrate --noinput 

# Collect all the static files
docker-compose exec web python manage.py collectstatic --no-input

```

## Start and stop BiGAnts-web
All the important files are stored inside persistant volumes and are available even after restarting the containers

**Start service**

Navigate to the repository folder
```shell script
docker-compose up -d
```

**Stop service**

Navigate to the repository folder

```shell script
docker-compose down
```

## Managing volumes and data
The docker containers store their data in three volumes which can be backup or moved when the application needs to 
move to another client.

**Backing up volumes**

**Deleting volumes / all data**

CAREFULL, docker does not ask for confirmation when deleting volumes!
```shell script
docker-compose down --volumes
```
