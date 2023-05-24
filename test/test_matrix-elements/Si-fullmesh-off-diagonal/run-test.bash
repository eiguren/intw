#!/usr/bin/env bash

# Run available scripts
for script in $(find ./ -name "matrix-element-test.bash");
do
  echo $(dirname "${script}")
  cd $(dirname "${script}")
  bash matrix-element-test.bash
  if [ $? -ne 0 ]; then
    exit 1
  fi
  cd ..
done
