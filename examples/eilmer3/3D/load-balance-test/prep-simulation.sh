#!/bin/bash

e3prep.py --job=test
e3loadbalance.py --job=test -n 16

