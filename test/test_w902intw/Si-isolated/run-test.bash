#!/usr/bin/env bash

# Run available scripts
for script in $(find ./ -name "w902intw-isolated-test.bash");
do
  echo $(dirname "${script}")
  cd $(dirname "${script}")
  bash w902intw-isolated-test.bash
  if [ $? -ne 0 ]; then
    exit 1
  fi
  cd ..
done
