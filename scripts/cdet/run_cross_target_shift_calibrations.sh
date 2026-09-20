#!/bin/bash
set -u

# Fit shift_ns only for the newly configured ECal-singles cross-target runs.
# Established reference runs 4344, 5710, and 5992 are intentionally omitted.
# Runs 5976 and 5979 are intentionally excluded because their timing is
# inconsistent with the ECal-singles trigger model.

runs=(
  3573 3575
  3602 3603 3604 3648 3682 3685 3686 3694
  3757 3758 3786 3788 3844 3845 3846 4003 4006
  4345 4400 4722 4735 4926 4929 4985 5066 5068
  5290 5292 5294 5423
  5722 5723 5724 5726 5794 5795
)

log_dir="CDet_cross_target_shift_logs"
backup_dir="CDet_cross_target_shift_preexisting_20260920"
failed_dir="CDet_cross_target_shift_failed"
summary="CDet_cross_target_shift_calibrations.tsv"

mkdir -p "$log_dir" "$backup_dir" "$failed_dir"
printf 'run\tstatus\tshift_ns\tconfig\trun_file\tlog\n' > "$summary"

if [[ ! -r CDet_calibration_dt.dat ]]; then
  echo "ERROR: CDet_calibration_dt.dat is missing." >&2
  exit 2
fi

overall_status=0
for run in "${runs[@]}"; do
  config="CDet_run${run}_projection.conf"
  run_file="CDet_run${run}.dat"
  log="$log_dir/CDet_run${run}_shift.log"
  backup="$backup_dir/$run_file"
  had_run_file=0

  if [[ ! -r "$config" ]]; then
    echo "Run $run: missing $config" | tee "$log"
    printf '%s\tFAILED_MISSING_CONFIG\t\t%s\t%s\t%s\n' "$run" "$config" "$run_file" "$log" >> "$summary"
    overall_status=1
    continue
  fi

  if [[ -e "$run_file" ]]; then
    had_run_file=1
    if [[ -e "$backup" ]]; then
      echo "Run $run: refusing to overwrite existing backup $backup" | tee "$log"
      printf '%s\tFAILED_EXISTING_BACKUP\t\t%s\t%s\t%s\n' "$run" "$config" "$run_file" "$log" >> "$summary"
      overall_status=1
      continue
    fi
    cp -p "$run_file" "$backup"
  fi

  echo "Run $run: fitting shift_ns with $config"
  root_status=0
  root -l -b -q "Run_CDet_Calibrate_RunTimingShift.C+(${run},-1,30.0,5.0,\"${config}\",false)" > "$log" 2>&1 || root_status=$?

  if [[ $root_status -eq 0 ]] &&
     grep -q '^  closure Gaussian-core centroid:' "$log" &&
     ! grep -q '^\[Run timing-shift calibration\] ERROR:' "$log" &&
     [[ -r "$run_file" ]]; then
    shift_ns=$(awk '$1 == "shift_ns" { value=$2 } END { print value }' "$run_file")
    if [[ -n "$shift_ns" ]]; then
      printf '%s\tPASS\t%s\t%s\t%s\t%s\n' "$run" "$shift_ns" "$config" "$run_file" "$log" >> "$summary"
      echo "Run $run: PASS, shift_ns=$shift_ns"
      continue
    fi
  fi

  overall_status=1
  if [[ $had_run_file -eq 1 ]]; then
    cp -p "$backup" "$run_file"
  elif [[ -e "$run_file" ]]; then
    mv "$run_file" "$failed_dir/$run_file"
  fi
  printf '%s\tFAILED\t\t%s\t%s\t%s\n' "$run" "$config" "$run_file" "$log" >> "$summary"
  echo "Run $run: FAILED; inspect $log"
done

echo "Summary: $summary"
exit "$overall_status"
