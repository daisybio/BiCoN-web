from django.urls import path

from .views import IndexView, AnalysisSetupView, documentation, source

app_name = 'clustering'
urlpatterns = [
    path('', IndexView.as_view(), name='index'),
    path('analysis/', AnalysisSetupView.as_view(), name='analysis_setup'),
    path('analysis/submit', AnalysisSetupView.as_view(), name='submit_job'),
    path('documentation/', documentation, name='documentation'),
    path('source/', source, name='source'),
]