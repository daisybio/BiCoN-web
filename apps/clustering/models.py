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
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True, blank=True)
    session_id = models.CharField(max_length=33, blank=True)
    submit_time = models.DateTimeField(editable=False)
    status = models.CharField(
        max_length=9,
        choices=JOB_STATES,
    )
    finished_time = models.DateTimeField(null=True, blank=True)
    # ppi_json = models.TextField(blank=True)
    ppi_json = models.FileField(upload_to=f'clustering/ppi/')
    # survival_plotly = models.TextField(blank=True)
    survival_plotly = models.FileField(upload_to=f'clustering/survival/')
    ppi_png = models.FileField(upload_to=f'clustering/ppi/')
    heatmap_png = models.FileField(upload_to=f'clustering/heatmap/')
    convergence_png = models.FileField(upload_to=f'clustering/convergence/')

    def __str__(self):
        return f'Job ID: {self.job_id} submitted on {self.submit_time}'

