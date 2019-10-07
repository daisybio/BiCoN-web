from django.urls import path

from django.conf.urls import include, url
from . import views_old as views
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

urlpatterns = [
    url(r'infopage.html', views.infopage),
    url(r'sources.html', views.sources),
    url(r'errorpage.html', views.errorpage),
    url(r'clustering.html', views.clustering),
    url(r'clustering_step_1.html', views.clustering_step_1),
    url(r'clustering_step_2.html', views.clustering_step_2),
    url(r'clustering_no_sessions.html', views.clustering_no_sessions),
    url(r'login.html', views.login_2),
    url(r'logout.html', views.logout_2),
    url(r'signup.html', views.signup),
    url(r'delete_user.html', views.delete_user),
    url(r'^infopage/$', views.infopage, name='infopage'),
    url(r'^sources/$', views.sources, name='sources'),
    url(r'^errorpage/$', views.errorpage, name='errorpage'),
    url(r'^on_one_page/$', views.clustering, name='clustering_6_part_1_3'),
    url(r'^on_two_pages/$', views.clustering_step_1, name='clustering_6_part_1_2'),
    url(r'^results/$', views.clustering_step_2, name='clustering_6_part_1_1'),
    url(r'^clustering/$', views.clustering, name='clustering'),
    url(r'^clustering_step_1/$', views.clustering_step_1, name='clustering_step_1'),
    url(r'^clustering_step_2/$', views.clustering_step_2, name='clustering_step_2'),
    url(r'^clustering_no_sessions/$', views.clustering_no_sessions, name='clustering_no_sessions'),
    url(r'^clustering_6_part_4/$', views.clustering, name='clustering_6_4'),
    url(r'^clustering_6_part_3/$', views.clustering_step_2, name='clustering_6_part_3_2'),
    url(r'^clustering_6_part_1/$', views.clustering_step_1, name='clustering_6_part_1'),
    url(r'^clustering_6/$', views.clustering_no_sessions, name='clustering_6'),
]

app_name = 'clustering_old'
urlpatterns += staticfiles_urlpatterns()
