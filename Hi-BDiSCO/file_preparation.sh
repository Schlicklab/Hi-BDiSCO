#!/bin/bash


[ "$(ls -A output/ini_struct)" ] && rm -r output/ini_struct/*

python src/analysis.py

