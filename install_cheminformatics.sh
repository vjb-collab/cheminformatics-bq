#!/bin/bash

## Create New Dataset to store the cheminformatics functions

bq mk --description "Dataset that will contain the cheminformatics functions" --dataset "cheminformatics" 

## Create Connection

bq mk --connection --display_name="Cheminformatics Connection" --connection_type=CLOUD_RESOURCE --location=US "cheminformatics-connection"

SERVICE_ACCOUNT=$(bq show --location=US --format=prettyjson --connection "cheminformatics-connection" | jq -r '.cloudResource.serviceAccountId')

## Create Cloud functions (should be looping over all functions under cloud-functions and calling /*/install_cf.sh)

PERM="roles/cloudfunctions.invoker"

cd cloud-functions/rdkit-pattern-fingerprint/

ENDPOINT=$(gcloud functions deploy rdkit-pattern-fingerprint --entry-point rdkit_pattern_fingerprint --runtime python39 \
    --trigger-http --format=json | jq -r '.httpsTrigger.url')

#EP=$(gcloud functions describe rdkit-pattern-fingerprint --format=json | jq -r '.httpsTrigger.url')

gcloud functions add-iam-policy-binding "rdkit-pattern-fingerprint" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

#printf 'CREATE FUNCTION cheminformatics.rdkit_pattern_fingerprint(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = "%s")' "${EP}" 

bq query --use_legacy_sql=false 'CREATE FUNCTION cheminformatics.rdkit_pattern_fingerprint(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = "https://us-central1-life-sciences-333615.cloudfunctions.net/rdkit-pattern-fingerprint")'

cd ../../


