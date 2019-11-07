"""BiGAnts_web URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
# Django imports
from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import include, path
from django.conf import settings

urlpatterns = [
    path('', include('apps.clustering.urls')),
    path('admin/', admin.site.urls),
]

# Serve media files during debugging through the django server aswell
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

# TEMPLATE TODO check if can be used
# urlpatterns = [
#     # Examples:
#     # url(r'^blog/', include('blog.urls', namespace='blog')),
#
#     # provide the most basic login/logout functionality
#     url(r'^login/$', auth_views.login,
#         {'template_name': 'core/login.html'}, name='core_login'),
#     url(r'^logout/$', auth_views.logout, name='core_logout'),
#
#     # enable the admin interface
#     url(r'^admin/', admin.site.urls),
# ]
