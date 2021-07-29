from celery.states import FAILURE, IGNORED, PENDING, RECEIVED, RETRY, REVOKED, STARTED, SUCCESS
from django.contrib.auth.models import User
from django.db import models


class Job(models.Model):
    SUBMITTED = 'SUBMITTED'
    JOB_STATES = [
        (SUBMITTED, 'submitted'),
        (FAILURE, 'failed'),
        (IGNORED, 'ignored'),
        (PENDING, 'pending'),
        (RECEIVED, 'received'),
        (RETRY, 'retry'),
        (REVOKED, 'revoked'),
        (STARTED, 'started'),
        (SUCCESS, 'success')
    ]

    job_id = models.UUIDField(primary_key=True, editable=False)
    job_name = models.CharField(max_length=30, default='', blank=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True, blank=True)
    session_id = models.CharField(max_length=33, blank=True)
    submit_time = models.DateTimeField(editable=False)
    status = models.CharField(
        max_length=9,
        choices=JOB_STATES,
    )
    finished_time = models.DateTimeField(null=True, blank=True)
    # ppi_json = models.TextField(blank=True)
    result_csv = models.FileField(upload_to=f'clustering/result/', blank=True)
    ppi_json = models.FileField(upload_to=f'clustering/ppi/')
    # survival_plotly = models.TextField(blank=True)
    survival_plotly = models.FileField(upload_to=f'clustering/survival/')
    ppi_png = models.FileField(upload_to=f'clustering/ppi/')
    heatmap_png = models.FileField(upload_to=f'clustering/heatmap/')
    convergence_png = models.FileField(upload_to=f'clustering/convergence/')
    netex_json = models.FileField(upload_to=f'clustering/netex/', null=True)

    def __str__(self):
        return f'Job ID: {self.job_id} submitted on {self.submit_time}, finished on {self.finished_time}'


class PpiNetworkCache(models.Model):
    # Network id (same as NDEx UUID)
    network_id = models.CharField(editable=False, max_length=50, primary_key=True, serialize=False)
    # Date and time this cached version has been last modified/updated
    last_modified = models.DateTimeField(auto_now=True)
    # Date and time the NDEx PPI network has been modified
    data_last_modified = models.DateTimeField()
    # The actual cached network data as string
    network_string = models.TextField()

    def __str__(self):
        return f'Network: {self.network_id} last modified on {self.last_modified}. NDEx last modified on {self.data_last_modified}'
