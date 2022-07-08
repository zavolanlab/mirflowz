#!/bin/bash

echo -e "JOB ID\t$SLURM_JOBID"
echo "=============================="
echo -e "rule\t$JOB_NAME"
echo "=============================="

{exec_job}

echo -e "==============================\n"
