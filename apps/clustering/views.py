from django.views.generic import TemplateView
from django.http import HttpResponse


class IndexView(TemplateView):
    template_name = "clustering/index.html"


class AnalysisSetupView(TemplateView):
    template_name = "clustering/analysis_setup.html"


def documentation(request):
    return HttpResponse("SOME DOCUMENTATION")


def source(request):
    return HttpResponse("SOURCE")
