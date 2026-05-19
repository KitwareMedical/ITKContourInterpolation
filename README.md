# Migrated to ITK main

The `MorphologicalContourInterpolation` module has been ingested into the Insight Toolkit (ITK) main repository and is no longer maintained here.

- **In-tree location:** [`Modules/Filtering/MorphologicalContourInterpolation`](https://github.com/InsightSoftwareConsortium/ITK/tree/main/Modules/Filtering/MorphologicalContourInterpolation)
- **Ingest PR:** [InsightSoftwareConsortium/ITK#6209](https://github.com/InsightSoftwareConsortium/ITK/pull/6209)
- **Merge commit:** [`bf2088d6a3`](https://github.com/InsightSoftwareConsortium/ITK/commit/bf2088d6a3d41c908d0eea99f2ee05bb7d11a259)

Future development, bug fixes, and pull requests should target ITK main. This repository's history is preserved for `git blame` and archival reference.

**Maintainer action:** please mark this repository **Archived** under repo Settings -> Danger Zone after merging this PR. See the [ITK ingestion strategy](https://github.com/InsightSoftwareConsortium/ITK/blob/main/Utilities/Maintenance/RemoteModuleIngest/INGESTION_STRATEGY.md) for context.

---

## Original module description

# ITKMorphologicalContourInterpolation

![Build, test, package](https://github.com/KitwareMedical/ITKMorphologicalContourInterpolation/workflows/Build,%20test,%20package/badge.svg)

[![PyPI](https://img.shields.io/pypi/v/itk-morphologicalcontourinterpolation.svg)](https://pypi.python.org/pypi/itk-morphologicalcontourinterpolation)

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/KitwareMedical/ITKMorphologicalContourInterpolation/blob/master/LICENSE)

## Overview

An [ITK](https://itk.org)-based implementation of morphological contour interpolation based off the paper:

Albu AB1, Beugeling T, Laurendeau D.
"A morphology-based approach for interslice interpolation of anatomical slices from volumetric images."
IEEE Trans Biomed Eng.
2008 Aug;55(8):2022-38.
doi: 10.1109/TBME.2008.921158.

👨‍💻 [`Live API Demo`] ✨

Documentation can be found in the [Insight Journal article](https://www.insight-journal.org/browse/publication/977):

Zukić D., Vicory J., McCormick M., Wisse L., Gerig G., Yushkevich P., Aylward S.
"ND Morphological Contour Interpolation",
The Insight Journal. January-December, 2016.
https://hdl.handle.net/10380/3563
https://insight-journal.org/browse/publication/977

- 🕮 [`Wasm Python Documentation`] 📚
- 🕮 [`JavaScript Documentation`] 📚

## Installation

### Python

To install the wasm Python packages:

```sh
pip install itkwasm-morphological-contour-interpolation
```

To install the native Python packages:

```sh
python -m pip install --upgrade pip
python -m pip install itk-morphologicalcontourinterpolation
```

### JavaScript

Install the JavaScript/TypeScript package:

```sh
npm install @itk-wasm/morphological-contour-interpolation
```

### C++

Since ITK 4.11.0, this module is available in the ITK source tree as a remote module. To enable it, set:

```
Module_MorphologicalContourInterpolation:BOOL=ON
```

in ITK's CMake build configuration.

## License

This software is distributed under the Apache 2.0 license. Please see the *LICENSE* file for details.

## Acknowledgements

This work is supported by NIH grant R01 EB014346, "Continued development and maintenance of the ITK-SNAP 3D image segmentation software."

[`Live API Demo`]: https://kitwaremedical.github.io/ITKContourInterpolation/ts/app/
[`Wasm Python Documentation`]: https://kitwaremedical.github.io/ITKContourInterpolation/py/docs/
[`JavaScript Documentation`]: https://kitwaremedical.github.io/ITKContourInterpolation/ts/docs/
