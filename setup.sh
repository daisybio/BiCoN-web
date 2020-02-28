#!/usr/bin/env bash

# ======================================================================================
# title			:  setup.sh
# description	:  Setting the .env file up for docker-compose
# usage			:  bash setup.sh [--configure] [--migrate-db] [--collect-static] [--update]
# author       	:  Kevin Yuan
# date			:  25.01.2020
# version		:  1.0
# notes			:  <Some notes about this script, dependencies>
# ======================================================================================

# Terminate this bash script as soon as an error occures
set -eu

# ================================== MISC Codebits ====================================
# Setting some colors
black=$(tput setaf 0)
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
redBg=$(tput setab 1)
greenBg=$(tput setab 2)
yellowBg=$(tput setab 3)
whiteBg=$(tput setab 7)
reset=$(tput sgr0)

termwidth="$(($(tput cols) / 2))"

# Centered printing with Padding
function center() {
  padding="$(printf '%0.1s' ={1..500})"
  printf '%*.*s %s %*.*s\n' 0 "$(((termwidth - 2 - ${#1}) / 2))" "$padding" "$1" 0 "$(((termwidth - 1 - ${#1}) / 2))" "$padding"
}

# Checking if flag (command line option) is set
# Taken from here https://stackoverflow.com/questions/2875424/correct-way-to-check-for-a-command-line-flag-in-bash
has_param() {
  local term="$1"
  shift
  for arg; do
    if [[ $arg == "$term" ]]; then
      return 0
    fi
  done
  return 1
}

wrap_text() {
  echo -e "$1" | fold -w $termwidth -s
}

# ================================= Argument Parser ===================================

# Welcome message
center ""
center "BiCoN web setup script"
center ""
echo -e

# Argument flags
run_all=true
configure=false
migrate_db=false
collect_static=false

# Run configuration file generation
if has_param '--configure' "$@"; then
  run_all=false
  configure=true
fi

# Migrate changes in the new database
if has_param '--migrate-db' "$@"; then
  run_all=false
  migrate_db=true
fi

# Collect new static files
if has_param '--collect-static' "$@"; then
  run_all=false
  collect_static=true
fi

# Update == Collect new static files && Migrate changes to the database
if has_param '--update' "$@"; then
  run_all=false
  migrate_db=true
  collect_static=true
fi

# If no option has been found, run everything and print options
if [[ $run_all == "true" ]]; then
  configure=true
  migrate_db=true
  collect_static=true

  # Print selection
  wrap_text "${yellow}No flag has been set. Running everything (configuring the containers; migrating the database; collecting all the static files).${reset}"
else
  echo -e "${yellow}The following operations will be perfomed"
  echo -e "- Configuring the containers from BiCoN-web.conf (creating the .env file): ${configure}"
  echo -e "- Migrating the database container based on the supplied migrations from django: ${migrate_db}"
  echo -e "- Collecting all the static files from the templates into the static directory: ${collect_static}"
  echo -e "${reset}"
fi

echo -e


# =========================== Configuring the container =============================

if [ $configure = "true" ]; then

  # Print section
  center "Configuring the container"
  echo -e

  # Check if config file exists
  echo -e "- Load BiCoN-web.conf"
  if ! test -f "BiCoN-web.conf"; then
    echo -e "${red}BiCoN-web.conf does not exist."
    echo -e "Exiting now...${reset}"
    exit 1
  fi
  # Load the config file (carefull, code in this file will be executed!)
  source BiCoN-web.conf

  # Backup old .env file if exists
  echo -e "- Check if configuration file already exists"
  if test -f ".env"; then
    echo -e "  + Old .env configuration file found, rename it into '.env.bak'"
    mv .env .env.bak
  fi

  echo -e
  echo -e "- Creating the new '.env' file:"
  : >.env

  # ------------ Django & Celery
  echo -e "--- Django & Celery:"
  echo -e "# --- Django & Celery ---" >>.env

  # Check for production
  echo -e "     + Django production: \c"

  if [[ $DJANGO_PRODUCTION == "true" ]]; then
    echo "DJANGO_SETTINGS_MODULE=BiCoN_web.settings.production" >>.env
    echo -e "${green}true (in production mode)${reset}"
  elif [[ $DJANGO_PRODUCTION == "false" ]]; then
    echo "DJANGO_SETTINGS_MODULE=BiCoN_web.settings.development" >>.env
    echo -e "${yellow}false (in debug mode)${reset}"
  else
    echo -e "${red}Error! Unknown parameter \"${DJANGO_PRODUCTION}\" - Valid options are \"true\" or \"false\""
    echo -e
    echo -e "Exiting now ...${reset}"
    exit 1
  fi

  # Allowed hosts
  echo -e "     + Django allowed hosts: \c"

  if [[ "$DJANGO_ALLOWED_HOSTS" =~ ^\[(\"|\').*(\"|\')\]$ ]]; then
    echo "DJANGO_ALLOWED_HOSTS=${DJANGO_ALLOWED_HOSTS}" >>.env
    echo -e "${green}${DJANGO_ALLOWED_HOSTS}${reset}"
  else
    echo -e "${red}Error! Your input: ${DJANGO_ALLOWED_HOSTS} is not valid. Please check if your input is wrapped in quotes"
    echo -e
    echo -e "Exiting now ...${reset}"
    exit 1
  fi

  # ------------ Postgres
  echo -e "--- Postgres:"
  echo -e "\n# --- Postgres ---" >>.env

  # Username
  echo -e "     + Postgres username: \c"

  if [ -z $POSTGRES_USER ]; then
    echo -e "${yellow}postgres (variable not set - default username)${reset}"
    POSTGRES_USER=postgres
  else
    echo -e "${green}${POSTGRES_USER}${reset}"
  fi

  echo "POSTGRES_USER=${POSTGRES_USER}" >>.env

  # Password
  echo -e "     + Postgres password: \c"

  if [ -z $POSTGRES_PASSWORD ]; then
    echo -e "${yellow}postgres (variable not set - default password)${reset}"
    POSTGRES_PASSWORD=postgres

  else
    echo -e "${green}*****${reset}"
  fi

  echo "POSTGRES_PASSWORD=${POSTGRES_PASSWORD}" >>.env

  # Database
  echo -e "     + Postgres database: \c"

  if [ -z $POSTGRES_DB ]; then
    echo -e "${yellow}postgres (variable not set - default database)${reset}"
    POSTGRES_DB=postgres
  else
    echo -e "${green}${POSTGRES_DB}${reset}"
  fi

  echo "POSTGRES_DB=${POSTGRES_DB}" >>.env

  # ------------ RabbitMQ
  echo -e "--- RabbitMQ:"
  echo -e "\n# --- RabbitMQ ---" >>.env

  # Username
  echo -e "     + RabbitMQ username: \c"
  if [ -z $RABBITMQ_DEFAULT_USER ]; then
    echo -e "${yellow}admin (variable not set - default username)${reset}"
    RABBITMQ_DEFAULT_USER=admin
  else
    echo -e "${green}${POSTGRES_USER}${reset}"
  fi

  echo "RABBITMQ_DEFAULT_USER=${RABBITMQ_DEFAULT_USER}" >>.env

  # Password
  echo -e "     + RabbitMQ password: \c"
  if [ -z $RABBITMQ_DEFAULT_PASS ]; then
    echo -e "${yellow}mypass (variable not set - default password)${reset}"
    RABBITMQ_DEFAULT_PASS=mypass
  else
    echo -e "${green}*****${reset}"
  fi
  echo "RABBITMQ_DEFAULT_PASS=${RABBITMQ_DEFAULT_PASS}" >>.env

  # ------------ NGINX
  echo -e "--- Nginx:"
  echo -e "\n# --- NGINX ---" >>.env

  # Published path
  echo -e "     + Nginx published path: \c"
  if [ -z $NGINX_PUBLISHED_PATH ]; then
    echo -e "${yellow}/ (variable not set - default path)${reset}"
  else
    # Remove leading or trailing slashes
    NGINX_PUBLISHED_PATH=${NGINX_PUBLISHED_PATH#/}
    NGINX_PUBLISHED_PATH=/${NGINX_PUBLISHED_PATH%/}

    echo -e "${green}${NGINX_PUBLISHED_PATH}${reset}"
  fi

  echo "NGINX_PUBLISHED_PATH=${NGINX_PUBLISHED_PATH}" >>.env

  # Port
  echo -e "     + Nginx published port: \c"
  if [ -z $NGINX_PUBLISHED_PORT ]; then
    echo -e "${yellow}8081 (variable not set - default port)${reset}"
    NGINX_PUBLISHED_PORT=8081
  else
    echo -e "${green}${NGINX_PUBLISHED_PORT}${reset}"
  fi

  echo "NGINX_PUBLISHED_PORT=${NGINX_PUBLISHED_PORT}" >>.env

  # --- END
  echo -e "\n${green}Generating .env from config file successfull${reset}"

fi

# =========================== Starting the containers =============================
# Container need to start / been build when we want to reconfigure them
if [[ $migrate_db == "true" ]] || [[ $collect_static == "true" ]]; then
  center "Build and start the docker containers"
  docker-compose up -d --build
fi

# =========================== Migrate the database =============================
if [[ $migrate_db == "true" ]]; then
  center "Migrate the database"
  docker-compose exec web python manage.py migrate --noinput
fi

# =========================== Collect all the static files =============================
if [[ $collect_static == "true" ]]; then
  center "Collect all the static files"
  docker-compose exec web python manage.py collectstatic --no-input
fi
