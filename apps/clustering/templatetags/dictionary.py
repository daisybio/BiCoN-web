from django import template

register = template.Library()


@register.simple_tag
def add_key(d, key, value):
    d[key] = value
    return d
