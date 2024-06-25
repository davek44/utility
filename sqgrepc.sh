#!/bin/bash

scancel `squeue --format "%.10i %.9P %.20j %.8u %.2t %.10M %.6D %Q %R" -u drk | grep $1 | awk '{print $1}' | xargs`
