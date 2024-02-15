# @itk-wasm/morphological-contour-interpolation

- 👨‍💻 [`Live API Demo`](https://kitwaremedical.github.io/ITKMorphologicalContourInterpolation/ts/app/) ✨
- 🕮 [`JavaScript Documentation`](https://kitwaremedical.github.io/ITKMorphologicalContourInterpolation/ts/docs/) 📚

[![npm version](https://badge.fury.io/js/@itk-wasm%2Fmorphological-contour-interpolation.svg)](https://www.npmjs.com/package/@itk-wasm/morphological-contour-interpolation)

> Morphology-based approach for interslice interpolation of anatomical slices from volumetric images.

## Installation

```sh
npm install @itk-wasm/morphological-contour-interpolation
```

## Usage

### Browser interface

Import:

```js
import {
  morphologicalContourInterpolation,
  setPipelinesBaseUrl,
  getPipelinesBaseUrl,
} from "@itk-wasm/morphological-contour-interpolation"
```

#### morphologicalContourInterpolation

*Interpolates contours between slices.*

```ts
async function morphologicalContourInterpolation(
  inputImage: Image,
  options: MorphologicalContourInterpolationOptions = {}
) : Promise<MorphologicalContourInterpolationResult>
```

|   Parameter  |   Type  | Description     |
| :----------: | :-----: | :-------------- |
| `inputImage` | *Image* | The input image |

**`MorphologicalContourInterpolationOptions` interface:**

|           Property          |             Type            | Description                                                                                                                                                                                        |
| :-------------------------: | :-------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|           `label`           |           *number*          | The label to interpolate. Interpolates all labels if set to 0 (default).                                                                                                                           |
|            `axis`           |           *number*          | Interpolate only along this axis. Interpolates along all axes if set to -1 (default).                                                                                                              |
|    `noHeuristicAlignment`   |          *boolean*          | Heuristic alignment of regions for interpolation is faster than optimal alignment.                                                                                                                 |
|   `noUseDistanceTransform`  |          *boolean*          | Using distance transform instead of repeated dilations to calculate the median contour is slightly faster, but produces lower quality interpolations.                                              |
|  `useCustomSlicePositions`  |          *boolean*          | Use custom slice positions (not slice auto-detection).                                                                                                                                             |
|     `noUseExtrapolation`    |          *boolean*          | Perform extrapolation for branch extremities. Branch extremities are defined as regions having no overlap with any region in the next slice. Extrapolation helps generate smooth surface closings. |
| `useBallStructuringElement` |          *boolean*          | Use ball instead of default cross structuring element for repeated dilations.                                                                                                                      |
|  `labeledSliceIndicesAxis`  |           *number*          | Axis along which the labeled slice indices are defined. Default is -1 (that is, auto-detection).                                                                                                   |
|  `labeledSliceIndicesLabel` |           *number*          | Label of the slice indices. Default is 1.                                                                                                                                                          |
|    `labeledSliceIndices`    |          *number[]*         | List of labeled slice indices. Default is empty.                                                                                                                                                   |
|         `webWorker`         | *null or Worker or boolean* | WebWorker for computation. Set to null to create a new worker. Or, pass an existing worker. Or, set to `false` to run in the current thread / worker.                                              |
|           `noCopy`          |          *boolean*          | When SharedArrayBuffer's are not available, do not copy inputs.                                                                                                                                    |

**`MorphologicalContourInterpolationResult` interface:**

|    Property   |   Type   | Description                     |
| :-----------: | :------: | :------------------------------ |
| `outputImage` |  *Image* | The output image                |
|  `webWorker`  | *Worker* | WebWorker used for computation. |

#### setPipelinesBaseUrl

*Set base URL for WebAssembly assets when vendored.*

```ts
function setPipelinesBaseUrl(
  baseUrl: string | URL
) : void
```

#### getPipelinesBaseUrl

*Get base URL for WebAssembly assets when vendored.*

```ts
function getPipelinesBaseUrl() : string | URL
```


### Node interface

Import:

```js
import {
  morphologicalContourInterpolationNode,
} from "@itk-wasm/morphological-contour-interpolation"
```

#### morphologicalContourInterpolationNode

*Interpolates contours between slices.*

```ts
async function morphologicalContourInterpolationNode(
  inputImage: Image,
  options: MorphologicalContourInterpolationNodeOptions = {}
) : Promise<MorphologicalContourInterpolationNodeResult>
```

|   Parameter  |   Type  | Description     |
| :----------: | :-----: | :-------------- |
| `inputImage` | *Image* | The input image |

**`MorphologicalContourInterpolationNodeOptions` interface:**

|           Property          |    Type    | Description                                                                                                                                                                                        |
| :-------------------------: | :--------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|           `label`           |  *number*  | The label to interpolate. Interpolates all labels if set to 0 (default).                                                                                                                           |
|            `axis`           |  *number*  | Interpolate only along this axis. Interpolates along all axes if set to -1 (default).                                                                                                              |
|    `noHeuristicAlignment`   |  *boolean* | Heuristic alignment of regions for interpolation is faster than optimal alignment.                                                                                                                 |
|   `noUseDistanceTransform`  |  *boolean* | Using distance transform instead of repeated dilations to calculate the median contour is slightly faster, but produces lower quality interpolations.                                              |
|  `useCustomSlicePositions`  |  *boolean* | Use custom slice positions (not slice auto-detection).                                                                                                                                             |
|     `noUseExtrapolation`    |  *boolean* | Perform extrapolation for branch extremities. Branch extremities are defined as regions having no overlap with any region in the next slice. Extrapolation helps generate smooth surface closings. |
| `useBallStructuringElement` |  *boolean* | Use ball instead of default cross structuring element for repeated dilations.                                                                                                                      |
|  `labeledSliceIndicesAxis`  |  *number*  | Axis along which the labeled slice indices are defined. Default is -1 (that is, auto-detection).                                                                                                   |
|  `labeledSliceIndicesLabel` |  *number*  | Label of the slice indices. Default is 1.                                                                                                                                                          |
|    `labeledSliceIndices`    | *number[]* | List of labeled slice indices. Default is empty.                                                                                                                                                   |

**`MorphologicalContourInterpolationNodeResult` interface:**

|    Property   |   Type  | Description      |
| :-----------: | :-----: | :--------------- |
| `outputImage` | *Image* | The output image |
