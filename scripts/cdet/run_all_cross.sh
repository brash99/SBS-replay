#!/bin/bash
N_GROUPS=$1

root -l -b -q "run_all_cross_runs.C($N_GROUPS)"