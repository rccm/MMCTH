#!/bin/bash
python ./tests/test_process_main.py &
PID=$!

MAX_MEM=0
while kill -0 $PID 2>/dev/null; do
    MEM=$(ps -o rss= -p $PID)
    if [ $MEM -gt $MAX_MEM ]; then
        MAX_MEM=$MEM
    fi
    sleep 0.1
done

echo "Maximum memory used: $((MAX_MEM / 1024)) MB"