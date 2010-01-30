#!/bin/bash

epydoc -v -u http://code.google.com/p/bayesian-inference --parse-only --html --no-frames BIP
cd BIP/Docs
make html latex
cd build/latex
make all-pdf
cd ../../../..
