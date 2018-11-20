#!/bin/sh
CXX=clang++ \
meson out \
-Db_sanitize=address -Db_lundef=false
