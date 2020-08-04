#!/bin/bash

# Author: Dogacan S. Ozturk

# Python installation script for all python packages required
# to run the HIME Framework.

# Python dependencies
easy_install pip
pip install --upgrade ipython
pip install --upgrade jupyter
pip install --upgrade numpy
pip install --upgrade matplotlib
pip install --upgrade scipy
pip install spacepy
pip install apexpy
pip install tables
