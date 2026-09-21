#!/usr/bin/env bash
# bin/resumable_download.sh URL OUT [MAX_ATTEMPTS]
set -uo pipefail
url="$1"; out="$2"; max="${3:-20}"

for attempt in $(seq 1 "$max"); do
    if wget -c -T 60 -O "$out" "$url"; then
        # wget reports complete; verify integrity
        if gzip -t "$out"; then
            exit 0
        fi
        echo "Attempt $attempt: download complete but corrupt; restarting from scratch" >&2
        rm -f "$out"
    else
        size=$(stat -c %s "$out" 2>/dev/null || echo 0)
        echo "Attempt $attempt/$max interrupted at $size bytes; resuming in 30s" >&2
    fi
    sleep 30
done

echo "ERROR: download failed after $max attempts: $url" >&2
exit 1
