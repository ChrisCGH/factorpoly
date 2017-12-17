#!/bin/bash

polynomial="$1"

curl -s http://localhost:8080/function/factorpoly -X POST --data-binary "${polynomial}" | firefox -new-tab "data:text/html;base64,$(base64)"
