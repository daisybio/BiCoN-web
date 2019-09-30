# BiGAnts Web - Code Notes

Some code analysis notes and page ideas / proposals

## URLS

+ `	url(r'infopage.html', views.infopage),`

  Infopage with guide. Set to default. "Home" in header

+ `	url(r'sources.html', views.sources),`

  Sources for application. "Sources" in Header

+ `	url(r'errorpage.html', views.errorpage),`

  Errorpage for redirection when errors occur.

+ ` url(r'on_one_page.html', views.clustering),`

  Default clustering page where everything is in one page (livereload?!)

+ `url(r'on_two_pages.html', views.clustering_step_1), `

  Analysis run (submission form) (but after redirect same as in views.clustering)

+ `	url(r'results.html', views.clustering_step_2),`

  Result page

+ `url(r'clustering_no_sessions.html', views.clustering_no_sessions),`

  In theory the same as in clustering, all in one page

## Layout

### Template

Create a master template for all pages.

Includes:

+ Main CSS / JS elements

+ Header
+ Footer

### Header

Header for every page. 

| Home     | Analysis                  | Documentation                   | Source                  | Login                   |
| -------- | ------------------------- | ------------------------------- | ----------------------- | ----------------------- |
| Infopage | Clustering  Analysis Page | Guide on how to use BiGAnts web | Sources for application | Login page with Sign up |



### Infopage

+ What is BiGAnts (short abstract?)
+ Links to other pages 

### Documentation

- How to run the algorithm (maybe add an overview with links?)
- Example Run + Parameters

### Analysis

Default page should be the one-page analysis.

+ Add short description
+ Add redirect to different analysis types (multi-page analysis & no session analysis)
+ Run Analysis Form (maybe link to the documentation in each step?)

Maybe default should be "multi page" but with automatic redirects? This way implementing a loading screen would be easier? And page is cleaner

### Sources

Sources for app and algorithm

### Login

Login Page for users.

Also provide registration form for new user sign up