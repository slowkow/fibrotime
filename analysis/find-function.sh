#!/usr/bin/env bash

rg -w "$1" -g '!website/*' -l # -g '!R/002_helper*'
