#!/bin/sh
find . -name \*.h -o -name \*.cpp -print | xargs grep TODO
