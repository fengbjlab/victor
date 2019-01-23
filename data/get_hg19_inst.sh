#! /usr/bin/env bash

ggID='1zhgbpU6EManvG1ZfCyIY45bOvejLNf27'
ggURL='https://drive.google.com/uc?export=download'
filename="$(curl -sc /tmp/gcokie "${ggURL}&id=${ggID}" | grep -o '="uc-name.*</span>' | sed 's/.*">//;s/<.a> .*//')"
getcode="$(awk '/_warning_/ {print $NF}' /tmp/gcokie)"
curl -Lb /tmp/gcokie -C - --insecure "${ggURL}&confirm=${getcode}&id=${ggID}" | tar xz
