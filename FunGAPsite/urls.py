"""FunGAPsite URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
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
from django.contrib import admin
from django.conf.urls import url, include
from fungap import views as fungap_views

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^homepage/', fungap_views.homepage),
    url(r'^captcha/', include('captcha.urls')),
    url(r'^confirm/', fungap_views.user_confirm),
    url(r'^new/', fungap_views.new),
    url(r'^custom/', fungap_views.custom),
    url(r'^delete_custom/', fungap_views.delete_custom),
    url(r'^delete_project/', fungap_views.delete_project),
    url(r'^get_markdown/', fungap_views.get_markdown),
    url(r'^get_log/', fungap_views.get_log),
]
