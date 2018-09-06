#!/bin/bash

while true
do
	memory_used=$( free | awk '/Mem/ {print $3}' )
	memory_used_gb=$(($memory_used / 1000000 ))
	echo $memory_used_gb
	sleep 1m
done
