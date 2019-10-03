# 1. Getting Started

You can run the application by running `docker-compose`, which starts the server and all necessities, by writing in a terminal (in the same directory as this file):

```
sudo allowed_hosts=['*'] cel_url='amqp://admin:mypass@rabbit:5672' secret_key='9z5(_w$5&=_)eve^u(--xcg%ge3dxi38m^d$yqol5#*atybvt6' docker-compose up
```

instead of `['*']`, which allows all hosts, you can write a number of the port that django is using, for example `['0.0.0.0']`. It is important to specify (even if it is *, meaning all addresses) allowed hosts since otherwise Django will not be able to show the website to you. The cel_url parameter (and port address) should be the same on all computers but it is possible to specify password and username (which will be described in 3. How the code works).

You can find out which port Django is using by running

```bash
sudo allowed_hosts=['*'] cel_url='amqp://admin:mypass@rabbit:5672'  secret_key='9z5(_w$5&=_)eve^u(--xcg%ge3dxi38m^d$yqol5#*atybvt6' docker-compose up
```

and waiting until a line appears in the terminal that says

> clust_app     | Starting development server at http://0.0.0.0:8000/

In this case, the address where the server runs is `0.0.0.0`. This is the address you need later to access the site.

IN ADDITION, in order to allow the live updates of algorithm progress, open a terminal in the same directory as this file (not on the Docker container, but on your local computer) and type:

```bash
python manage.py livereload 
```

OR, to use your own secret key (my_secret_key):

```
export SECRET_KEY="my_secret_key"
python manage.py livereload
```

in the same terminal.

This gets the server running. You can then access the website by typing:

> http://0.0.0.0:8000/clustering/clustering.html

And substitute `0.0.0.0:8000` by the address and port number that is printed out by the application as above.

# 2. How to use the web application

To run a bi-clustering analysis on the web application, you can use either built-in datasets or your own data. You must upload or select data by clicking through a multi-step form.
To make an example run, click the big button that says "run new analysis" and then check the checkbox that says "use pre-compiled data". Select one of the possibilites in the dropdown menu (it selects expression datasets).
On the next step, click "import ndex" and select an option from the dropdown menu (it contains PPI networks for the analysis). Skip the third step. On the fourth step, type "10" in the first field and "20" in the second field. This selects a target gene number of 10 to 20 for the resulting networks. Click submit and the analysis is run. You will see a progress image that tells you the current progress of your analysis. If it says "done", scroll down to view the tabs with results. You can select which part of the result to show by selecting from the big clickable bar.

Own data must consist of one file with expression data in CSV or TSV format (the exact format is visible on the website), a file with PPI consisting of two columns with interaction partners, and optionally of clinical patient data in CSV format. If you want to analyze survival of patients, you must indicate which column of your clinical data file lists the survival time of patients.

You can change the exact settings for the algorithm in the fourth tab.

To run an enrichment analysis, scroll down to the results and select the tab "run enrichment". 



# 3. How the code works 

### What runs where:

There are 4 docker containers that each run a task or a server that is necessary for the web application:

- **Livereload** (named *clust_app_livereload_1*): A python-based application/library for Django needed for a good display of the loading page. Livereload makes every file reload automatically when it is changed.
There appears to still be a bug with livereload on Docker, so it is possible that this part of the server will be replaced by javascript code with the same functionality.
- **RabbitMQ** (*clust_app_rabbit_1*): Needed as server where celery executes its tasks. Must run in the background for celery to run properly.
- **Celery** (*celery_container*): Used as queuing system to execute tasks synchronously. Executes all computationally intensive tasks such as running the algorithm or processing results.
- **Django** (*clust_app*): The server for running and processing requests. This container is where the actual server runs and where paths to files can be referenced.

as well as a fifth container for the database needed for Django, called *clust_app_db_1*.

There is a shared volume that can be accessed by all containers under `/code`. This is needed to transfer files from one to another container by the server.

In order to show the progress of the algorithm while running, LiveReload is needed to run in both the container and the folder where the container is started. For this, open two tabs in the terminal, in the first run: 

```
(optionally, if you use a different secret key): export SECRET_KEY='YOUR_OWN_SECRET_KEY'
```

and

```
python manage.py livereload
```

In the second:

```
sudo docker compose-up
```

This ensures that everytime a gif or png changes, which is needed to display the running status of the algorithm, it is reloaded in the page displayed to the user.

In addition to the other processes in containers, a crontab runs on the clust_app_web container at /etc/crontab to automatically delete old unused user files (files with results, e.g heatmap png). It is run every 59 minutes, and removes old files from the userfiles-directory where result files are stored. The time after which files are deleted can be set via modifying the line 
```
*/59 * * * * find /clust_app/clustering/static/userfiles/* -mmin +59 -type f -delete
```
to

```
*/5 * * * * find /clust_app/clustering/static/userfiles/* -mmin +YOUR_TIME_IN_MINUTES -type f -delete
```

It can be extended to other directories by adding lines with a link to that directory, to the file.

If you want to change the password and username for RabbitMQ for Celery, you can substitute "admin" and "mypass" by your username and password in the command

```bash
sudo allowed_hosts=['*'] cel_url='amqp://admin:mypass@rabbit:5672'  secret_key='9z5(_w$5&=_)eve^u(--xcg%ge3dxi38m^d$yqol5#*atybvt6' docker-compose up
```



### Environment variables

In `docker-compose.yml`, general specifications are listed for the different parts of the application. There you can specify a different password for RabbitMQ by modifying the line RABBITMQ_DEFAULT_USER.
You can change the port on which the website is running by modifying the command

```
python manage.py runserver 0.0.0.0:8000 
```

and substituting 0.0.0.0 by another address.

Most of the settings for the server are in "clust_app/clust_app/settings.py", such as the communication between Django, Celery and RabbitMQ. There are comments which line is used to do what, and what must be modified to change a specific setting. It is possible to add default values for the input variables (e.g. ['*'] for accepted_hosts) in "clust_app/clust_app/settings.py" by changing

```
ALLOWED_HOSTS = os.getenv("DJANGO_ALLOWED_HOSTS")
```

to

```
ALLOWED_HOSTS = os.getenv("DJANGO_ALLOWED_HOSTS", ['*'])
```

Then, the input variables become optional, and the server will run with the default values if none are given.

### What is where (in the code):

Templates (HTML pages visible to the user) are connected in a file at `clustering/urls.py` with a function that processes the request for each template. If you want to change the link for one of the pages, add a line below "urlpatterns" in this form:

```
url(r'new_address.html', views.your_view),
```

with your_view as the name of the method in "clustering/views.py" that is used to process a request submitted via "new_address.html" (e.g. the method for clustering.html is called clustering).

Files with results are stored in a static directory for which the location is defined in "clust_app/clust_app/settings.py". This static directory is located in a shared volume (path is "/code") that is accessible for all containers (celery, livereload, django, rabbitmq). This is important because the celery container manipulates data that must be visible to the django container. In some cases, files are copied from the shared volume to the container running the Django server. 

The code needed for running a request is in three different files, in the folder "clustering":
- views.py - processes the request, will be described below.
- tasks.py - running the algorithm on the queuing system and pre/postprocessing files
- weighted_aco_lib - parts of the algorithm that are connected to a method in tasks.py

Additionally, for logged in users, two functions in "clustering/models.py" are used:
- save_user_data_3 - to save the uploaded files
- list_user_data - to list the stored input files
- list_user_data_2 - to list the stored result files

The files with results are stored in the following directories:
/code/clustering/static/userfiles - here (on the shared volume) they are stored first after being generated by Celery tasks
and
/clust_app/clustering/static/userfiles/ - on the clust_app_web container, where they can be referenced for the templates.
For every run of the algorithm, files are created with the heatmap, json data for PPI graph, html version of the survival plot, convergence plot, metadata (for later loading) and gene lists (for later enrichment analyses).



### How a request runs on the server:

If you submit a request, several steps are executed (on the clust_app container) in views.py. 
For the loading page, the gif with loading symbol is copied to its place so that it is visible in the page. The current status of the run is written and read from a log file.
Data are preprocessed by the file preprocessing methods in tasks.py so that CSV files are converted to TSV. It is checked whether a column in the expression data defines two clinical clusters; if this is the case, these clusters are later displayed in the heatmap. Automatically, a column containing information on pre-defined clinical clusters is recognized and assigned a standard title ("disease_type"). Patient identifiers are taken from the first column, and gene IDs from the first row. The PPI data are filtered for rows with two columns with IDs of protein interaction partners.
Two scripts are used to run the algorithm for the data in tasks.py: algo_output_task and script_output_task. These tasks are executed via Celery, and their results are written on the shared volume (under /code). In algo_output_task, the algorithm itself is run. A preprocessing script (aco_preprocessing_strings) in weighted_aco_lib.py is run, the script also automatically recognizes whether the expression data are log2-values (if they contain negative values) and ant runs are executed via the ants() method in weighted_aco_lib.py. The output data (clusters of patients and genes) are taken as an input in script_output_task and there are processed to output format (json and png). They are stored by the methods run in the celery container in a directory on the shared volume and copied by views.py, which is run on the django container, to the local static files directory.
If a logged-in user wants to display old results, there is a list of stored result files (by date). The filename of the expression file is generated based on the current date when saving, and is written as the value of a dropdown option in the input form when selecting. 

The version of the page with separate pages for input and results ("on_two_pages.html") redirects the user to the result page after running the analysis. A parameter in the cache is set to "done" so that if the user hits reload or runs an enrichment analysis, he stays on the result page. The user can go back to the input page by clicking "start new analysis". If the user presses this button, an form is submitted with an input parameter "new analysis", and in the method in views.py, the "done" parameter is set to false and the user is redirected to the input page.



### How sessions work:

For user managment, Session IDs (stored in the cache) are used. When a user first accesses the page, he is assigned a new session ID and that ID is stored in the cache. All files that are used for progress display and all files with results contain the session ID in their path (e.g. progress_[SESSION_ID].png).
When a constructor method is called, the session ID is read from the cache and included in the file paths. It is passed to the methods that write these files; if no session IDs are used, the parameter for session ID is set to "none" which is recognized by these methods. The rendered template contains links to all files as variables.
Updates to these files are detected automatically by the livereload application and reflected to the user.



### How the cache works:

In order to keep result data visible to the user after hitting reload or running an enrichment analysis, the metadata (if existing) and paths to the result files are kept as variables in the cache. If a user submits an "empty" request (e.g. by hitting the reload button), or runs an enrichment analysis, the views.py script will look up whether these variables are in the cache, and returns a page with a link to the corresponding result files. The Session ID is also stored as a variable in the cache.
Currently, the cache timeout is set as 10000 seconds to enable long sessions for users. This can be changed (e.g. to five minutes) in "settings.py" by modifying the line below "CACHES=" to (e.g.) "timeout: 300".



### How the templates work:

The results are referenced in the templates by setting django variables (in the views.py script) with paths to json files for PPI data and images for the other plots. In a js function, a Sigma.js instance is generated and loads the PPI network on the path referenced by the django variable. The survival plot is referenced as a html file that is set as data source for a html-"object" tag. Metadata and enrichment results are read as a table from dictionaries passed by the Django function in "views.py" and heatmap and convergence plot are referenced as static png files.
Two javascript libraries are used for the templates (they are included in the Docker container), sigma.js for displaying PPI graphs and dataTable.js for a good display of enrichment results in tables.



### Which methods do what:

In views.py, methods are defined that run code after the user accesses a web page.

There are three different constructors that render a web page for the clustering application with slightly different functionality:

- clustering_no_sessions: A static version of the webpage that does not use session IDs and with input and results on the same page.
- clustering_step_1 (accessible under "http://0.0.0.0:8000/clustering/on_two_pages.html": A version of the page with different pages for input and result display. Redirects automatically to the result page after completion of the algorithm. Uses Session IDs (stored in the cache) for user management.
- clustering (accessible under ("http://0.0.0.0:8000/clustering/on_one_page.html"): Version of the page with input and results on the same page. Uses Session IDs (in the cache) for user management.

# 4. Troubleshooting:
If you want to run bash commands on one of the containers (e.g. for troubleshooting), you can select the ID of the container by typing

```
sudo docker ps
```

and you get the ID of the container (YOUR_CONTAINER_ID) based on name. Then type 

```
sudo docker exec -it YOUR_CONTAINER_ID bash
```

and you can now run bash commands on the container.

Some Errors that might occur:

- Celery does not start properly and RabbitMQ produces an error upon starting:

After usable space on the server was totally full, this error can occur. It is possibly solved by typing in a shell on the clust_app_rabbit_1 container and running

```
sudo rabbitmqctl stop
sudo rabbitmq-server
```

If this does not help, then a log file of the RabbitMQ server was corrupted. This is no big deal and can be solved by running (on a shell in the clust_app_rabbit_1 container)

```
mkdir -p /tmp/badrabbit/ 
sudo mv /var/lib/rabbitmq/mnesia/rabbit[[[COMPUTER_NAME]]]/queues/* /tmp/badrabbit/;
```

substitute COMPUTER_NAME by the name of your computer/the server in the shell, that is shown in the terminal at the start of each line.

substitute COMPUTER_NAME by the name of your computer/the server in the shell, that is shown in the terminal at the start of each line.

- RabbitMQ does not do anything after running a request:

If this happens, look for a line in the output of RabbitMQ that says that the memory is insufficient. RabbitMQ needs at least 50 MB of free memory on the system, otherwise it will not take any requests from Celery. Your request will be run again as soon as there is more than 50 MB of free memory.

- The database gives a "duplicate table" error on the shell:
e.g. "psycopg2.errors.DuplicateTable: relation "clustering_graphform" already exists"
This is a database error that is produced because there are two containers running postgres. This error can happen if you accidentally built the postgres contatiner twice. Delete the older postgres
container to avoid this error.



