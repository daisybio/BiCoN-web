# Some codesnippets

## Fix Paths
Replace absolute media paths like `"/code/clustering/static(.*)"`
```
"/code/clustering/static([^"])*"
```
to
```
path.join(settings.MEDIA_ROOT, 'clustering/userfiles$1')
```


## Button Grou

```html
        <div class="btn-group">
            <a href="#about" class="btn btn-outline-secondary btn-lg">Read more</a>
            <a href="{% url 'clustering:analysis_setup' %}" class="btn btn-primary btn-lg">Start analysis!</a>
        </div>
```