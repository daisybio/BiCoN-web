 

# BiGAnts Web - ToDo

## Code cleanup

+ [ ] Reduce Pages (e.g. `clustering_6_part_4`) in `urls.py`

  + [ ] One overview page of "clustering", default, one page. Advance options for no seassion and two pages?
  + [ ] Set clustering infopage as default
  + [ ] Remove .html suffix and remove clustering subpage

+ [ ] Standardize Layout of header (Home Analaysis Guide?!) 

  Home(views.infopage), Analysis(views.clustering), Sources(views.sources), Sign up

+ [ ] Check redirects (?) `^clustering_6_part_1/$`

## Bug Fixes

+ [x] Analysis does not run: Userfiles could not be found

  Create folder `clustering/static/userfiles/`

  Update Docker?

+ [ ] Enrichment does not work

  => Works in `clustering_step_1.html`

  Breaks in clustering, why?

+ [ ] User  sign up does not work

  ```
  clust_app     | [03/Sep/2019 22:21:06] "GET /clustering/signup.html%20target= HTTP/1.1" 200 4462
  clust_app     | Internal Server Error: /clustering/signup.html target=
  clust_app     | Traceback (most recent call last):
  clust_app     |   File "/usr/local/lib/python3.6/site-packages/django/core/handlers/exception.py", line 34, in inner
  clust_app     |     response = get_response(request)
  clust_app     |   File "/usr/local/lib/python3.6/site-packages/django/core/handlers/base.py", line 126, in _get_response
  clust_app     |     response = self.process_exception_by_middleware(e, request)
  clust_app     |   File "/usr/local/lib/python3.6/site-packages/django/core/handlers/base.py", line 124, in _get_response
  clust_app     |     response = wrapped_callback(request, *callback_args, **callback_kwargs)
  clust_app     |   File "/clust_app/clustering/views.py", line 81, in signup
  clust_app     |     os.mkdir(userdir)
  clust_app     | FileNotFoundError: [Errno 2] No such file or directory: 'user_uploaded_files/test'
  clust_app     | [03/Sep/2019 22:21:17] "POST /clustering/signup.html%20target= HTTP/1.1" 500 44288
  ```

+ [ ] No session run does not work

  ```
  celery_container | Traceback (most recent call last):
  celery_container |   File "/usr/local/lib/python3.6/site-packages/billiard/process.py", line 327, in _bootstrap
  celery_container |     self.run()
  celery_container |   File "/usr/local/lib/python3.6/site-packages/billiard/process.py", line 114, in run
  celery_container |     self._target(*self._args, **self._kwargs)
  celery_container |   File "/clust_app/clustering/weighted_aco_lib.py", line 133, in ant_job_paral
  celery_container |     tot_score,gene_groups,patients_groups,new_scores,wars,no_int = ant_job(GE,N,H,th,clusters,probs,a,b,cost,m,n,patients,count_big,cost_limit,L_g_min,L_g_max,G,ge,seed)
  
  ```

+ [ ] Fix opening "result page" before running analysis ? Or if user is not signed in



____

##### Less critical

+ [ ] Secure Postgres DB 

  ```
  db_1          | WARNING: No password has been set for the database.
  db_1          |          This will allow anyone with access to the
  db_1          |          Postgres port to access your database. In
  db_1          |          Docker's default configuration, this is
  db_1          |          effectively any other container on the same
  db_1          |          system.
  db_1          | 
  db_1          |          Use "-e POSTGRES_PASSWORD=password" to set
  db_1          |          it in "docker run".
  ```

+ [ ] Secure celary & Django?

  ```
  celery_container | [2019-09-03 21:57:50,523: INFO/MainProcess] mingle: all alone
  celery_container | [2019-09-03 21:57:50,601: WARNING/MainProcess] /usr/local/lib/python3.6/site-packages/celery/fixups/django.py:202: UserWarning:
  celery_container | 
  celery_container | Using settings.DEBUG leads to a memory leak, never use this setting in production environments!
  celery_container | [2019-09-03 21:57:50,601: INFO/MainProcess] celery@095af29d04bb ready.
  ```



## Additional features

+ [ ] Only submit same job once
+ [ ] Progress bar
+ [ ] Better layout (fixed elements somehow?)



## Helpful commands:

### Connect empty Docker container to folder

**CAREFUL RUNNING AS ROOT ON HOST SYSTEM!!!!!!**

```bash
docker run -v /data/home/users/kyuan/projects/BigAnts:/root/BigAnts -i -t --rm ubuntu bash
```

Will automatically be removed after usage. Uses an empty Ubuntu base container.

