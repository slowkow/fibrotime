#!/usr/bin/env bash

des="ig:/var/www/immunogenomics.io/public/fibrotime"

#rsync -rlzuvh \
rsync -ahvz \
   --delete \
   --no-perms --no-owner --no-group --ignore-times --omit-dir-times \
   --progress \
   --exclude=*.swp \
   --exclude=archive \
   --exclude=data/*.R \
   --exclude=.git \
   --exclude=launch.sh \
   --exclude=*.Rproj \
   --exclude=.* \
   "./" "$des"
