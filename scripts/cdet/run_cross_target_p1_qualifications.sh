#!/bin/bash
set -u

# Frozen-sample p1 qualification for every accepted ECal-singles cross-target
# run. QUALIFIED=YES and QUALIFIED=NO are both valid scientific outcomes.
# Runs 5711 and 6077 are LH2; runs 5976 and 5979 have an HCal-carried trigger
# condition. All four are intentionally excluded.

runs=(
  3573 3575
  3602 3603 3604 3648 3682 3685 3686 3694
  3757 3758 3786 3788 3844 3845 3846 4003 4006
  4344 4345 4400 4722 4735 4926 4929 4985 5066 5068
  5290 5292 5294 5423
  5710 5722 5723 5724 5726 5794 5795
  5992
)

log_dir="CDet_cross_target_p1_qualification_logs"
artifact_backup_dir="CDet_cross_target_p1_prior_artifacts_20260920"
run_file_backup_dir="CDet_cross_target_p1_run_file_backups_20260920"
summary_tsv="CDet_cross_target_p1_qualifications.tsv"

mkdir -p "$log_dir" "$artifact_backup_dir" "$run_file_backup_dir"
printf '%s\n' 'run	status	residual_p1	residual_p1_error	significance	candidate_p1	timing_effect_ns	resolution_ratio	gate_discovery	gate_layers	gate_practical	gate_closure	gate_folds	gate_resolution	config	summary	log' > "$summary_tsv"

if [[ ! -r CDet_calibration_dt.dat ]]; then
  echo "ERROR: CDet_calibration_dt.dat is missing." >&2
  exit 2
fi

read_field() {
  awk -v key="$2" '$1 == key { print $2; exit }' "$1"
}

technical_failures=0
for run in "${runs[@]}"; do
  config="CDet_run${run}_projection.conf"
  run_file="CDet_run${run}.dat"
  sample_file="CDet_run${run}_frozen_p1_samples.tsv"
  result_file="CDet_run${run}_frozen_p1_qualification.txt"
  log="$log_dir/CDet_run${run}_p1_qualification.log"
  run_file_backup="$run_file_backup_dir/$run_file"

  if [[ ! -r "$config" || ! -r "$run_file" ]]; then
    echo "Run $run: missing configuration or timing file" | tee "$log"
    printf '%s\tERROR_MISSING_INPUTS\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n' "$run" "$config" "$result_file" "$log" >> "$summary_tsv"
    technical_failures=1
    continue
  fi

  if [[ -e "$run_file_backup" ]]; then
    echo "Run $run: refusing to overwrite $run_file_backup" | tee "$log"
      printf '%s\tERROR_EXISTING_BACKUP\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n' "$run" "$config" "$result_file" "$log" >> "$summary_tsv"
    technical_failures=1
    continue
  fi
  cp -p "$run_file" "$run_file_backup"
  before_checksum=$(shasum -a 256 "$run_file" | awk '{print $1}')

  for artifact in "$sample_file" "$result_file"; do
    if [[ -e "$artifact" ]]; then
      mv "$artifact" "$artifact_backup_dir/$artifact"
    fi
  done

  echo "Run $run: frozen-sample p1 qualification"
  root_status=0
  root -l -b -q "Run_CDet_Qualify_RunSpecificP1_FrozenSample.C+(${run},-1,\"${config}\",5,3.0,0.5,2.0,0.5,1.05)" > "$log" 2>&1 || root_status=$?

  after_checksum=$(shasum -a 256 "$run_file" | awk '{print $1}')
  if [[ "$after_checksum" != "$before_checksum" ]]; then
    cp -p "$run_file_backup" "$run_file"
    echo "Run $run: ERROR; qualification modified $run_file; original restored" >> "$log"
    printf '%s\tERROR_RUN_FILE_MODIFIED\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n' "$run" "$config" "$result_file" "$log" >> "$summary_tsv"
    technical_failures=1
    continue
  fi

  qualified=""
  if [[ -r "$result_file" ]]; then
    qualified=$(read_field "$result_file" qualified)
  fi
  if [[ $root_status -ne 0 || ( "$qualified" != "0" && "$qualified" != "1" ) ]] ||
     grep -q '^\[Frozen p1\] ERROR:' "$log"; then
    printf '%s\tERROR\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n' "$run" "$config" "$result_file" "$log" >> "$summary_tsv"
    echo "Run $run: ERROR; inspect $log"
    technical_failures=1
    continue
  fi

  if [[ "$qualified" == "1" ]]; then
    status="QUALIFIED"
  else
    status="NOT_QUALIFIED"
  fi
  residual=$(read_field "$result_file" residual_p1)
  residual_error=$(read_field "$result_file" residual_p1_error)
  significance=$(read_field "$result_file" full_significance)
  candidate=$(read_field "$result_file" candidate_p1)
  timing_effect=$(read_field "$result_file" timing_effect_ns)
  resolution_ratio=$(read_field "$result_file" resolution_ratio)
  gate_discovery=$(read_field "$result_file" gate_discovery)
  gate_layers=$(read_field "$result_file" gate_layers)
  gate_practical=$(read_field "$result_file" gate_practical)
  gate_closure=$(read_field "$result_file" gate_closure)
  gate_folds=$(read_field "$result_file" gate_folds)
  gate_resolution=$(read_field "$result_file" gate_resolution)

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$run" "$status" "$residual" "$residual_error" "$significance" "$candidate" "$timing_effect" "$resolution_ratio" "$gate_discovery" "$gate_layers" "$gate_practical" "$gate_closure" "$gate_folds" "$gate_resolution" "$config" "$result_file" "$log" >> "$summary_tsv"
  echo "Run $run: $status (residual=$residual +/- $residual_error)"
done

echo "Summary: $summary_tsv"
exit "$technical_failures"
