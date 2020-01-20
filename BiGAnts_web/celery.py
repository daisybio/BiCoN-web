import os

from celery import Celery

# TODO: Add other file? or method to change settings.dev
# set the default Django settings module for the 'celery' program.
os.environ.get("DJANGO_SETTINGS_MODULE", "BiGAnts_web.settings.development")

app = Celery('BiGAnts-web')

# Using a string here means the worker doesn't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')

# Load task modules from all registered Django app configs.
app.autodiscover_tasks()


# @app.task(bind=True)
# def debug_task(self):
#     print('Request: {0!r}'.format(self.request))