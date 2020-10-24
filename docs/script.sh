#!/bin/bash
sphinx-autogen -o reference/api ./reference/*.rst; make clean; make html