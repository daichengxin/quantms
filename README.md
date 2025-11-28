# quantms - ReadTheDocs Branch

This branch contains only the documentation configuration for ReadTheDocs.

## About quantms

**bigbio/quantms** is a bioinformatics pipeline for Quantitative Mass Spectrometry (MS). For the full pipeline code and documentation, please visit:

- **Main Repository**: [https://github.com/bigbio/quantms](https://github.com/bigbio/quantms)
- **Documentation**: [https://docs.quantms.org](https://docs.quantms.org/en/latest/)
- **nf-core Page**: [https://nf-co.re/quantms](https://nf-co.re/quantms)

## Purpose of this Branch

This branch is used exclusively for ReadTheDocs configuration. It redirects visitors from `quantms.readthedocs.io` to the main documentation site at `docs.quantms.org`.

## Documentation Structure

```
docs/
├── source/
│   ├── conf.py      # Sphinx configuration with redirect settings
│   └── index.rst    # Landing page with redirect notice
├── requirements.txt # Python dependencies for building docs
└── README.md        # Documentation overview
```

## Building Documentation Locally

To build the documentation locally:

```bash
pip install -r docs/requirements.txt
cd docs
sphinx-build -b html source _build/html
```

## ReadTheDocs Configuration

The `.readthedocs.yaml` file configures the ReadTheDocs build environment:
- Uses Python 3.11
- Builds with Sphinx
- Installs dependencies from `docs/requirements.txt`

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
