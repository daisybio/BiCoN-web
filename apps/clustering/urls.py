from django.urls import path

from .views import IndexView, AnalysisSetupView, submit_analysis, DocumentationView, SourceView

app_name = 'clustering'
urlpatterns = [
    path('', IndexView.as_view(), name='index'),
    path('analysis/', AnalysisSetupView.as_view(), name='analysis_setup'),
    path('analysis/submit', submit_analysis, name='submit_job'),
    path('documentation/', DocumentationView.as_view(), name='documentation'),
    path('source/', SourceView.as_view(), name='source'),
]