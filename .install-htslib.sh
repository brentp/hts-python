#!/bin/bash
DIR=${TRAVIS_BUILD_DIR}
BRANCH=develop
git clone https://github.com/samtools/htslib.git ${DIR}/htslib
cd ${DIR}/htslib
git checkout $BRANCH
sudo make install
sudo ldconfig
