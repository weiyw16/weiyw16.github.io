#!/bin/bash

cp -r /Users/weiyw/Library/Containers/com.coderforart.MWeb3/Data/Documents/themes/Site/English/* ./
git add ./*
git commit -m $1
git push

