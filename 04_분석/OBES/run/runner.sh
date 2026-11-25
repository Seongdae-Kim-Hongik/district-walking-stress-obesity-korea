#!/bin/bash
# Job runner. Usage: RSLOT=<path to 01_스크립트/rslot.sh> MAXJ=4 runner.sh <log folder> <phase1.jobs> [<phase2.jobs> ...]
# One line per job: name<TAB>time limit in minutes<TAB>command. Jobs in a phase run in parallel (at most MAXJ, started when the 1-minute load is below 14, 20 s apart);
# phases run in order. Each job's exit code is written to rc_ledger.tsv; a job that exceeds its limit is stopped (exit code 124).
set -m
LOGD=$1; shift; mkdir -p "$LOGD"; LED="$LOGD/rc_ledger.tsv"; MAXJ=${MAXJ:-4}; RSLOT=${RSLOT:?RSLOT 경로 필요}
[ -f "$LED" ] || printf "phase\tjob\tstart\tend\telapsed_s\trc\tload1_at_start\tgovernor_pids\ttimeout_min\n" > "$LED"
et2s(){ local e=$1 d=0 h=0 m=0 s=0; e=${e// /}; if [[ $e == *-* ]]; then d=${e%%-*}; e=${e#*-}; fi
  IFS=: read -r a b c <<< "$e"; if [ -n "$c" ]; then h=$a; m=$b; s=$c; else m=$a; s=$b; fi
  echo $(( 10#$d*86400 + 10#$h*3600 + 10#$m*60 + 10#$s )); }
load1(){ sysctl -n vm.loadavg | awk '{print $2}'; }
for JF in "$@"; do
  PH=$(basename "$JF" .jobs); NAMES=(); CMDS=(); TOS=(); while IFS=$'\t' read -r nm to cmd; do [ -z "$nm" ] && continue; [[ $nm == \#* ]] && continue; NAMES+=("$nm"); TOS+=("$to"); CMDS+=("$cmd"); done < "$JF"
  echo "[runner $(date '+%F %T')] phase $PH: ${#NAMES[@]} jobs" >> "$LOGD/runner.log"
  PIDS=(); STARTS=(); next=0; ndone=0; declare -a DONE=()
  while [ $ndone -lt ${#NAMES[@]} ]; do
    for i in "${!PIDS[@]}"; do P=${PIDS[$i]}; [ -z "$P" ] && continue; [ "${DONE[$i]}" = 1 ] && continue
      if ! kill -0 "$P" 2>/dev/null; then wait "$P"; rc=$?; DONE[$i]=1; ndone=$((ndone+1)); en=$(date +%s)
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$PH" "${NAMES[$i]}" "$(date -r ${STARTS[$i]} '+%F %T')" "$(date -r $en '+%F %T')" $((en-STARTS[$i])) "$rc" "${LD[$i]}" "${GV[$i]}" "${TOS[$i]}" >> "$LED"
        echo "[runner $(date '+%F %T')] done ${NAMES[$i]} rc=$rc" >> "$LOGD/runner.log"
      else et=$(ps -o etime= -p "$P" 2>/dev/null); if [ -n "$et" ] && [ "$(et2s "$et")" -gt $(( ${TOS[$i]} * 60 )) ]; then
          echo "[runner $(date '+%F %T')] WATCHDOG kill ${NAMES[$i]} (etime $et > ${TOS[$i]} min)" >> "$LOGD/runner.log"; kill -TERM -- -"$P" 2>/dev/null; fi; fi
    done
    run=0; for i in "${!PIDS[@]}"; do [ -n "${PIDS[$i]}" ] && [ "${DONE[$i]}" != 1 ] && run=$((run+1)); done
    if [ $next -lt ${#NAMES[@]} ] && [ $run -lt $MAXJ ] && awk -v l="$(load1)" 'BEGIN{exit !(l<14)}'; then
      i=$next; LD[$i]=$(load1); GV[$i]=none
      ( cd "$(dirname "$RSLOT")" && SK_RSLOTS=$MAXJ "$RSLOT" bash -c "${CMDS[$i]}" ) > "$LOGD/${NAMES[$i]}.log" 2>&1 &
      PIDS[$i]=$!; STARTS[$i]=$(date +%s); next=$((next+1))
      echo "[runner $(date '+%F %T')] start ${NAMES[$i]} pid ${PIDS[$i]} load ${LD[$i]} governor ${GV[$i]}" >> "$LOGD/runner.log"; sleep 20; continue
    fi
    sleep 15
  done
  unset DONE
done
echo "[runner $(date '+%F %T')] ALL PHASES DONE" >> "$LOGD/runner.log"; echo RUNNER_DONE
