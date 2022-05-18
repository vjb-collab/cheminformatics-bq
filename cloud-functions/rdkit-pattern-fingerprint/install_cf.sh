#!/bin/bash

cd "$(dirname "$0")"

gcloud functions deploy rdkit-pattern-fingerprint --entry-point rdkit_pattern_fingerprint --runtime python39 --trigger-http --allow-unauthenticated