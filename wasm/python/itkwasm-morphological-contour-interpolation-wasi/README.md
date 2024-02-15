# itkwasm-morphological-contour-interpolation-wasi

[![PyPI version](https://badge.fury.io/py/itkwasm-morphological-contour-interpolation-wasi.svg)](https://badge.fury.io/py/itkwasm-morphological-contour-interpolation-wasi)

Morphology-based approach for interslice interpolation of anatomical slices from volumetric images. WASI implementation.

This package provides the WASI WebAssembly implementation. It is usually not called directly. Please use [`itkwasm-morphological-contour-interpolation`](https://pypi.org/project/itkwasm-morphological-contour-interpolation/) instead.


## Installation

```sh
pip install itkwasm-morphological-contour-interpolation-wasi
```

## Development

```sh
pip install pytest itkwasm-image-io
pip install -e .
pytest

# or
pip install hatch
hatch run test
```
