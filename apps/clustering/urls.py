from django.urls import path

from .views import IndexView, AnalysisSetupView, submit_analysis, analysis_status, analysis_result, results,\
    DocumentationView, SourceView, test

app_name = 'clustering'
urlpatterns = [
    path('', IndexView.as_view(), name='index'),
    # Analysis setup page (overview page / first page)
    path('analysis/', AnalysisSetupView.as_view(), name='analysis_setup'),
    # Only for post-requests. If called via head, just redirect to status
    path('analysis/submit', submit_analysis, name='submit_job'),
    # Page to monitor current status of job
    path('analysis/status/<uuid:analysis_id>', analysis_status, name='analysis_status'),
    # Page for looking at one specific result
    path('analysis/result/<uuid:analysis_id>', analysis_result, name='analysis_result'),
    # Page for all the results
    path('analysis/results', results, name='analysis_results'),

    # TEST PAGE
    path('analysis/test', test, name='test'),

    # Page for documentation and source
    path('documentation/', DocumentationView.as_view(), name='documentation'),
    path('source/', SourceView.as_view(), name='source'),
]