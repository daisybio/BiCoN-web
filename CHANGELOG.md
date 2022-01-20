# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
### Changed
- Switched the Docker image from Python 3.6 (EOL) to Python 3.10
- Fixed the default settings file so `manage.py` can be called without setting environmental variables
### Removed
### Deploy Notes

## [1.1.0] - 2021-02-21
Based on BiCoN version 1.2.14

### Added
- Include App version and used BiCoN version in the footer
- Caching of the PPI networks
- PpiNetworkCache model and database table  
- Usage documentation on the setup.sh script

### Changed
- Using BiCoN 1.2.14 now
- Rewritten import_ndex function for (12x - 20x) speed up of downloading PPI networks
- Updated the citation section in the documentation.html section as well
- Updated the about section to include more links and information about used packages
- Optimized all imports, removed unused ones
- Reduced requirements.txt to a handful of actually used ones (switch from hard == to compatible ~=)
- Bumped all dependencies to the newest ones if possible (most notably we are using Django 3.0)
- Fixed style for Plotly (4.0) to match the other styles
- Explicitly set "X-Frame-Origin" to "sameorigin"
- Restrict CSRF token to this explicit site only (allows for compatibility between multiple different tools on the same server)

### Removed
- Old requirements.txt in trash folder
- Unused / commented code

### Deploy notes
Due to the updated database structure a "migrate" command is necessary

## [1.0.0] - 2021-02-20
Initial release of BiCoN-web based on BiCoN version 1.2.3
### Added
- New release tag
### Changed
- Fixed the import of NDEx PPI networks to new format
### Removed
