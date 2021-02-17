from django import template
from django.conf import settings

register = template.Library()


@register.simple_tag
def bicon_web_version():
    return settings.BICON_WEB_VERSION_NUMBER


@register.simple_tag
def bicon_version():
    return settings.BICON_PACKAGE_VERSION_NUMBER
