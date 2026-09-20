#!/bin/bash
set -u

runs=(
  3573 3575 3602 3603 3604 3648 3682 3685 3686 3694
  3757 3758 3786 3788 3844 3845 3846 4003 4006 4345
  4400 4722 4735 4926 4929 4985 5066 5068 5290 5292
  5294 5423 5722 5723 5724 5726 5794 5795 5976 5979
)

summary="CDet_cross_target_ECalTimingWindow_PROPOSALS.tsv"
printf 'run\tstatus\tlog\tpng\tpdf\troot\n' > "$summary"

for run in "${runs[@]}"; do
  log="CDet_run${run}_ECalTimingWindow_PROPOSAL.log"
  png="CDet_run${run}_ECalTimingWindow_PROPOSAL.png"
  pdf="CDet_run${run}_ECalTimingWindow_PROPOSAL.pdf"
  root_file="CDet_run${run}_ECalTimingWindow_PROPOSAL.root"
  if root -l -b -q "Run_CDet_Select_ECalTimingWindow.C+(${run})" > "$log" 2>&1 &&
     test -s "$png" && test -s "$pdf" && test -s "$root_file"; then
    status="generated"
  else
    status="failed"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$run" "$status" "$log" "$png" "$pdf" "$root_file" >> "$summary"
done

generated=$(awk -F '\t' 'NR > 1 && $2 == "generated" {count++} END {print count+0}' "$summary")
failed=$(awk -F '\t' 'NR > 1 && $2 == "failed" {count++} END {print count+0}' "$summary")
printf 'Generated %s proposals; %s failures. Index: %s\n' \
  "$generated" "$failed" "$summary"
test "$failed" -eq 0
