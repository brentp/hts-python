#!/bin/bash
DIR=${TRAVIS_BUILD_DIR}
BRANCH={v2.20.0}
git clone https://github.com/arq5x/bedtools2.git ${DIR}/bedtools2
cd ${DIR}/bedtools2
git checkout $BRANCH
make
