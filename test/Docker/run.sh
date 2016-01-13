#!/bin/sh

script_dir="`cd $(dirname $0); pwd`"

docker run \
  --rm \
  -v $script_dir/../..:/usr/src/ITKMorphologicalContourInterpolation \
    insighttoolkit/morphologicalcontourinterpolation-test \
      /usr/src/ITKMorphologicalContourInterpolation/test/Docker/test.sh
