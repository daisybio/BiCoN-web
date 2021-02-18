# BiCoN-web

Web tool of the [BiCoN](https://github.com/biomedbigdata/BiCoN) package for _network-constrained biclustering of patients and multi-omics data_.

A hosted version and demo of BiCoN-web can be found under: https://exbio.wzw.tum.de/bicon/

## Table of Contents

* [Installing BiCoN-web](#installing-bicon-web)
  + [Configuration of BiCoN-web](#configuration-of-bicon-web)
    - [Configuration options](#configuration-options-)
  + [Deploying (automatically)](#deploying--automatically-)
  + [Deploying (manually)](#deploying--manually-)
* [Updating and maintaining BiCoN-web (`setup.sh`)](#updating-and-maintaining-bicon-web---setupsh--)
* [Starting and stoping BiCoN-web](#starting-and-stoping-bicon-web)
* [Managing volumes and data](#managing-volumes-and-data)
* [Cite](#cite)
* [Contact](#contact)

## Installing BiCoN-web
BiCoN-web comes in an easy to use `docker-compose` format.

To deploy BiCoN-web in a production environment, please ensure you have a working docker (and docker-compose) environment running. 
You find a manual for install docker on your operating system [here](https://docs.docker.com/install/).

Once docker is installed and running, clone this repository using the  following command:
```shell script
# First clone this repository and change into the created directory
git clone https://github.com/biomedbigdata/BiCoN-web.git && cd BiCoN-web
```

### Configuration of BiCoN-web

The docker-compose file for BiCoN-web parses its configuration from the `.env` file. A `.env.sample` file is provided but it's not recommended to edit the file directly. We suggest you to specify your configuration inside the  `BiCoN-web.conf` file and generate the `.env` file with the `setup.sh` script.

```shell script
# First create a copy of the configuration sample
cp BiCoN-web.conf.sample BiCoN-web.conf
```

Open the file in a text editor and change the configuration to match your use case.

#### Configuration options:

+ `DJANGO_PRODUCTION`: set to either `true` or `false` 
+ `DJANGO_ALLOWED_HOSTS`: Only really relevant for production use. Please check the Django documentation for more information
+ `POSTGRES/RABBITMQ` passwords and username: Credentials for PostgreSQL and RabbitMQ Docker containers. They can be left with the default passwords as the containers are not reachable outside of the docker-compose setup
+ `NGINX_PUBLISHED_PATH`: If BiCoN-web should not be deployed on the root of the domain (e.g. example.com) but rather example.com/bicon, set the variable to `bicon`
+ `NGINX_PUBLISHED_PORT`: The port on which the nginx reverse proxy will bind to

### Deploying (automatically)

The `setup.sh` script generates the needed `.env` file and deploys the whole docker instance for you.

```shell script
# Just execute the script without any arguments
deploy.sh

# Enjoy your instance of BiCoN-web
```

Your instance should be available on http://localhost (or the port you have specified)

### Deploying (manually)

To install BiCoN-web (running it for the first time) follow these steps:
```shell script
# To generate a valid .env file from your BiCoN-web.conf configuration execute:
setup.sh --configure

# Deploy and build the containers
docker-compose up -d --build

# Apply migrations to the database (make sure the containers are up an running, else you will get an error)
docker-compose exec web python manage.py migrate --noinput 

# Collect all the static files
docker-compose exec web python manage.py collectstatic --no-input

# Enjoy your instance of BiCoN-web
```

Your instance should be available on http://localhost (or the port you have specified)

## Updating and maintaining BiCoN-web (`setup.sh`)

The `setup.sh` is also useful for deploying updated versions of BiCoN-web on existing setups.

The following arguments are available:

```shell script
bash setup.sh [--configure] [--migrate-db] [--collect-static] [--update]
```

`No arguments specified`:

+ equals: `bash setup.sh --configure --migrate-db --collect-static`
+ Useful for new installations

`--update`

+ equals: `bash setup.sh --configure --migrate-db --collect-static`
+ Useful for updating an existing instance (with existing `.env`) file

`--migrate-db`

+ Executes the Django migrate command inside the docker container

`--collect-static`

+ Executes the Django collect static command instide the docker container

**Force update cached NDEx PPI networks**
The cache is automatically invalidated 24h after generation or if the NDEx network has changed.
In case you want to force update the cache execute this command:
```bash script
docker-compose exec web python manage.py update_ppi_cache
```

## Starting and stoping BiCoN-web

All the important files are stored inside persistent volumes and are available even after restarting the containers

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
docker-compose down --volumes  # THIS COMMAND DELETES ALL THE STORED DATA
```

## Cite

BiCoN was developed by the [Big Data in BioMedicine group](https://github.com/biomedbigdata/BiCoN/blob/master/biomedical-big-data.de) and [Computational Systems Medicine group](https://compsysmed.de/) at [Chair of Experimental Bioinformatics](https://www.baumbachlab.net/).

If you use BiCoN in your research, we kindly ask you to cite the following manuscript: `Olga Lazareva, Stefan Canzar, Kevin Yuan, Jan Baumbach, David B Blumenthal, Paolo Tieri, Tim Kacprowski*, Markus List*, BiCoN: Network-constrained biclustering of patients and omics data, Bioinformatics, 2020;, btaa1076, https://doi.org/10.1093/bioinformatics/btaa1076`

** joint last author

## Contact

If you have difficulties using BiCoN-web, please open an issue at out [GitHub](https://github.com/biomedbigdata/BiCoN-web/issues/new) repository and/or tag @kevihiiin.ypo