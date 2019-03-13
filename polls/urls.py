from django.urls import path


from django.conf.urls import include, url
from . import views
from . import views_2
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
app_name = 'polls'
urlpatterns = [
	#path('', views.IndexView.as_view(), name='index'),
	#path('<int:pk>/', views.DetailView.as_view(), name='detail'),
	#path('<int:pk>/results/', views.ResultsView.as_view(), name='results'),
	#path('<int:question_id>/vote/', views.vote, name='vote'),
#	url(r'home.html', views.home),
	url(r'homepage.html', views_2.homepage),
	#url(r'predictions_3.html', views_2.prediction),
#	url(r'list.html', views.list_2),
	url(r'infopage.html', views.infopage),
	url(r'errorpage.html', views.errorpage),
	#url(r'clustering_7.html', views.clustering_7),
	url(r'clustering_6_part_3.html', views.clustering_6_part_3_2),
	url(r'clustering_6_part_2_2.html', views.clustering_6_part_2_2),
	url(r'clustering_6_part_2.html', views.clustering_6_part_2),
	url(r'clustering_6_part_1.html', views.clustering_6_part_1),
	url(r'clustering_6.html', views.clustering_6),
	#url(r'clustering_5.html', views.clustering_5),
	#url(r'clustering_4.html', views.clustering_4),
	#url(r'clustering_3.html', views.clustering_3),
	#url(r'clustering_2.html', views.clustering_2),
	#url(r'clustering.html', views.clustering),
	url(r'login.html', views.login_2),
	url(r'logout.html', views.logout_2),
#	url(r'other_test.html', views.list_3),
#	url(r'test3.html', views.list_2),
#	url(r'^list/$', views.list_2, name='list'),
	url(r'^infopage/$', views.infopage, name='infopage'),
	url(r'^errorpage/$', views.errorpage, name='errorpage'),
	#url(r'^clustering_7/$', views.clustering_7, name='clustering_7'),
	url(r'^clustering_6_part_3/$', views.clustering_6_part_3_2, name='clustering_6_part_3_2'),
	url(r'^clustering_6_part_2_2/$', views.clustering_6_part_2_2, name='clustering_6_part_2_2'),
	url(r'^clustering_6_part_2/$', views.clustering_6_part_2, name='clustering_6_part_2'),
	url(r'^clustering_6_part_1/$', views.clustering_6_part_1, name='clustering_6_part_1'),
	url(r'^clustering_6/$', views.clustering_6, name='clustering_6'),
	#url(r'^clustering_5/$', views.clustering_5, name='clustering_5'),
	#url(r'^clustering_4/$', views.clustering_4, name='clustering_4'),
	#url(r'^clustering_3/$', views.clustering_3, name='clustering_3'),
	#url(r'^clustering_2/$', views.clustering_2, name='clustering_2'),
	#url(r'^clustering/$', views.clustering, name='clustering'),
#	url(r'^other_test/$', views.list_3, name='other_test'),
        #url(r'^upload/$', views.list_2, name='imageupload')
] 
urlpatterns += staticfiles_urlpatterns()
