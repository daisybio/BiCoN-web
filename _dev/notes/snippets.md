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