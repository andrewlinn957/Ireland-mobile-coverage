# Ireland Mobile Signal Map

Source backup of the national Ireland mobile coverage map.

- Live Site: https://ireland-mobile-signal-map.andrewlinn.chatgpt.site
- Saved Site version: 2
- Source snapshot commit: `ca3c1a0b209b14fc63b5bca75e489e36e86bf98f`
- Static map interface and current coverage data: `dist/`
- National coverage calculation: `scripts/build-national-coverage.py`
- Model inputs: `model-inputs/`

The saved build script uses 500 m projected grid spacing (about 300 m on the ground at Irish latitudes) and a 1.5 m receiver height. This repository captures the saved Site source snapshot; separate run logs, checkpoints, or later experimental outputs are not included.
