FROM debian:8
MAINTAINER Insight Software Consortium <community@itk.org>

RUN apt-get update && apt-get install -y \
  build-essential \
  curl \
  cmake \
  git \
  libexpat1-dev \
  libhdf5-dev \
  libjpeg-dev \
  libpng12-dev \
  libpython3-dev \
  libtiff5-dev \
  python \
  ninja-build \
  wget \
  vim \
  zlib1g-dev

WORKDIR /usr/src

# 2016-01-12
RUN git clone https://github.com/Kitware/KWStyle.git && \
  cd KWStyle && \
  git checkout ddf8765ff13577572ec468be7b93deac77abed60 && \
  mkdir ../KWStyle-build && \
  cd ../KWStyle-build && \
  cmake \
    -G Ninja \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DCMAKE_BUILD_TYPE=Release \
      ../KWStyle && \
  ninja && \
  ninja install && \
  cd .. && rm -rf KWStyle*

# 2016-07-29 master
ENV ITK_GIT_TAG c7e7c79d0ebd2ee0dd69b34da3d0752a2e859dd8
RUN git clone git://itk.org/ITK.git && \
  cd ITK && \
  git checkout ${ITK_GIT_TAG} && \
  cd ../ && \
  mkdir ITK-build && \
  cd ITK-build && \
  cmake \
    -G Ninja \
    -DCMAKE_CXX_FLAGS='-std=c++11' \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON \
    -DITK_LEGACY_REMOVE:BOOL=ON \
    -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
    -DITK_USE_SYSTEM_LIBRARIES:BOOL=ON \
    -DModule_ITKBinaryMathematicalMorphology:BOOL=ON \
    -DModule_ITKTestKernel:BOOL=ON \
    -DModule_ITKDistanceMap:BOOL=ON \
    -DModule_ITKImageGrid:BOOL=ON \
    -DModule_ITKIOImageBase:BOOL=ON \
    ../ITK && \
  ninja && \
  find . -name '*.o' -delete

RUN git clone https://github.com/KitwareMedical/ITKRLEImage.git && \
  cd ITKRLEImage && \
  git checkout fa4a5aa5e37845b5b3d6cb6cbf69e75f2367f799 && \
  mkdir ../ITKRLEImage-build && \
  cd ../ITKRLEImage-build && \
  cmake \
    -G Ninja \
    -DCMAKE_CXX_FLAGS='-std=c++11' \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DCMAKE_BUILD_TYPE=Release \
    -DITK_DIR=/usr/src/ITK-build \
      ../ITKRLEImage && \
  ninja

RUN mkdir -p /usr/src/ITKMorphologicalContourInterpolation-build
