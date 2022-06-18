#!/bin/bash 

## Create New Dataset to store the cheminformatics functions

bq mk --description "Dataset that will contain the cheminformatics functions" --dataset "cheminformatics_test" 

## Create Connection

bq mk --connection --display_name="Cheminformatics Connection Test" --connection_type=CLOUD_RESOURCE --location=US "cheminformatics-connection-test"

## Get Service Account associated to Connection

SERVICE_ACCOUNT=$(bq show --location=US --format=prettyjson --connection "cheminformatics-connection-test" | jq -r '.cloudResource.serviceAccountId')

echo "Connection created with service account: ${SERVICE_ACCOUNT}"

## Give service account the cloud run invoker role (necessary for cloud functions gen2)

PROJ=$(gcloud config list --format 'value(core.project)')

gcloud projects add-iam-policy-binding $PROJ --quiet --member=serviceAccount:$SERVICE_ACCOUNT --role=roles/run.invoker

## Create Cloud functions 

PERM="roles/cloudfunctions.invoker"

## install rdkit-pattern-fingerprint-test

gcloud beta functions deploy rdkit-pattern-fingerprint-test --gen2 --region "us-east1" --entry-point rdkit_pattern_fingerprint_test --runtime python39 \
    --trigger-http --quiet --memory=512MB --timeout=240s --max-instances=9000  > /dev/null

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-pattern-fingerprint-test --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-pattern-fingerprint-test" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics_test.rdkit_pattern_fingerprint_test(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection-test` OPTIONS (endpoint = @url)'

## install rdkit-substructure-match

gcloud beta functions deploy rdkit-substructure-match --gen2 --region "us-east1" --entry-point rdkit_substructure_match --runtime python39 \
    --trigger-http --quiet --memory=512MB --timeout=240s --max-instances=9000  > /dev/null

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-substructure-match --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-substructure-match" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics_test.rdkit_substructure_match(fragment_smiles STRING, smiles STRING) RETURNS BOOL REMOTE WITH CONNECTION `us.cheminformatics-connection-test` OPTIONS (endpoint = @url)'

## install rdkit-molecular-descriptors

gcloud beta functions deploy rdkit-molecular-descriptors --gen2 --region "us-east1" --entry-point rdkit_molecular_descriptors --runtime python39 \
    --trigger-http --quiet --memory=512MB --timeout=240s --max-instances=9000  > /dev/null

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-molecular-descriptors --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-molecular-descriptors" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics_test.rdkit_molecular_descriptors(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection-test` OPTIONS (endpoint = @url)'

## install substructure search stored procedure

bq query --use_legacy_sql=false \
'CREATE OR REPLACE PROCEDURE `cheminformatics_test.substructure_search`(fragment_smiles STRING)
BEGIN
declare json_return STRING;
declare fragment_num_carbon, fragment_num_nitrogen, fragment_num_oxygen, fragment_num_fluorine, bit_count_fragment INT64;
declare fragment_pattern_fp BYTES;
declare fragment_bit_count INT64;
set json_return=`cheminformatics_test.rdkit_pattern_fingerprint_test`(fragment_smiles);
set fragment_pattern_fp = FROM_HEX(JSON_VALUE(json_return, "$.fp_pattern_long_as_binary_hex"));
set fragment_num_carbon = CAST(JSON_VALUE(json_return, "$.num_carbon") as INT64);
set fragment_num_oxygen = CAST(JSON_VALUE(json_return, "$.num_oxygen") as INT64);
set fragment_num_nitrogen = CAST(JSON_VALUE(json_return, "$.num_nitrogen") as INT64);
set fragment_num_fluorine = CAST(JSON_VALUE(json_return, "$.num_fluorine") as INT64);
set fragment_bit_count = bit_count(fragment_pattern_fp);
select smiles from `exports.savi_fp_clustered`
where
num_carbon >= fragment_num_carbon and
num_nitrogen >= fragment_num_nitrogen and
num_oxygen >= fragment_num_oxygen and
num_fluorine >= fragment_num_fluorine and
bit_count(fragment_pattern_fp & fp_pattern_as_binary) = fragment_bit_count
and
`cheminformatics_test.rdkit_substructure_match`(fragment_smiles, smiles);
END'

## install morgan similarity search 

bq query --use_legacy_sql=false \
'CREATE OR REPLACE PROCEDURE `cheminformatics_test.morgan_similarity`(smiles STRING)
BEGIN
declare json_return STRING;
declare morgan_fp BYTES;
declare in_num_carbon, in_num_nitrogen, in_num_oxygen, in_num_fluorine INT64;
set json_return=`cheminformatics_test.rdkit_pattern_fingerprint_test`(smiles);
set morgan_fp = FROM_HEX(JSON_VALUE(json_return, "$.fp_morgan_as_binary_hex"));
set in_num_carbon = CAST(JSON_VALUE(json_return, "$.num_carbon") as INT64);
set in_num_oxygen = CAST(JSON_VALUE(json_return, "$.num_oxygen") as INT64);
set in_num_nitrogen = CAST(JSON_VALUE(json_return, "$.num_nitrogen") as INT64);
set in_num_fluorine = CAST(JSON_VALUE(json_return, "$.num_fluorine") as INT64);
 select smiles from `exports.savi_fp_clustered`
where
num_carbon in (in_num_carbon-4, in_num_carbon - 3, in_num_carbon - 2, in_num_carbon - 1, in_num_carbon, in_num_carbon + 1, in_num_carbon + 2, in_num_carbon+3, in_num_carbon+4)
and
num_oxygen in (in_num_oxygen-4, in_num_oxygen - 3, in_num_oxygen - 2, in_num_oxygen - 1, in_num_oxygen, in_num_oxygen + 1, in_num_oxygen + 2, in_num_oxygen+3, in_num_oxygen+4)
and
bit_count(morgan_fp & fp_morgan_as_binary) / bit_count(morgan_fp || fp_morgan_as_binary) > 0.3;
END'

## wait one minute for permissions to propagate
echo "Waiting for permissions to propagate ..."
sleep 90

