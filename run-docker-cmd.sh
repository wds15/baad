#!/bin/bash

RUN=${1-bash}

echo "Running $RUN with stan_base-r4-baad:v2_19"

docker run --mount type=bind,source="$(pwd)",target=/app --workdir /app -it stan_base-r4-baad:v2_19 $RUN
