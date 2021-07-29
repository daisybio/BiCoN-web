# project imports
from .common import *

# include i18n internationalisation
from .i18n import *

# Dummy variable to import so PyCharm does not remove it while optimizing imports
__dummy = LANGUAGE_CODE

# ##### DEBUG CONFIGURATION ###############################
# turn off all debugging
DEBUG = True

# ##### DATABASE CONFIGURATION ############################
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': os.environ.get('POSTGRES_DB', 'postgres'),
        'USER': os.environ.get('POSTGRES_USER', 'postgres'),
        'PASSWORD': os.environ.get('POSTGRES_PASSWORD', 'postgres'),
        'HOST': 'db',
        'PORT': '5432',
    },
}

# ##### CACHE CONFIGURATION ############################
# CACHES = {
#     'default': {
#         'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
#         'LOCATION': 'my_cache_table',
#     }
# }

# ##### CELARY CONFIGURATION ############################
CELERY_BROKER_URL = f'amqp://{os.environ.get("RABBITMQ_DEFAULT_USER", "admin")}:{os.environ.get("RABBITMQ_DEFAULT_PASS", "mypass")}@rabbit:5672'
CELERY_TASK_SERIALIZER = 'pickle'
CELERY_RESULT_SERIALIZER = 'pickle'
CELERY_RESULT_BACKEND = "amqp"
CELERY_ACCEPT_CONTENT = ['pickle', 'json', 'msgpack', 'yaml']

# ##### APPLICATION CONFIGURATION #########################
INSTALLED_APPS = DEFAULT_APPS

# ##### SECURITY CONFIGURATION ############################
ALLOWED_HOSTS = ['*']

# TODO: Make sure, that sensitive information uses https
# TODO: Evaluate the following settings, before uncommenting them
# redirects all requests to https
# SECURE_SSL_REDIRECT = True
# session cookies will only be set, if https is used
# SESSION_COOKIE_SECURE = True
# how long is a session cookie valid?
# SESSION_COOKIE_AGE = 1209600

# validates passwords (very low security, but hey...)
# AUTH_PASSWORD_VALIDATORS = [
#    { 'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator', },
#    { 'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator', },
#    { 'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator', },
#    { 'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator', },
# ]

# the email address, these error notifications to admins come from
# SERVER_EMAIL = 'root@localhost'

# how many days a password reset should work. I'd say even one day is too long
# PASSWORD_RESET_TIMEOUT_DAYS = 1
