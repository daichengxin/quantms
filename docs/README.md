# bigbio/quantms: Documentation

This directory contains the ReadTheDocs configuration for the quantms project.

## About

The quantms documentation has moved to a new location. This ReadTheDocs configuration redirects visitors to the main documentation site.

- **Main Documentation**: [https://docs.quantms.org](https://docs.quantms.org/en/latest/)
- **nf-core Documentation**: [https://nf-co.re/quantms](https://nf-co.re/quantms)

## Structure

- `source/conf.py` - Sphinx configuration with redirect settings
- `source/index.rst` - Landing page with redirect notice
- `requirements.txt` - Python dependencies for building docs

## Building

To build locally:

```bash
pip install -r requirements.txt
sphinx-build -b html source _build/html
```
