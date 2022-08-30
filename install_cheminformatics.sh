#!/bin/bash 

## Create dataset to store the cheminformatics functions

bq mk  --force=true --description "Dataset that will contain the cheminformatics functions" --label=package:cheminformatics  --dataset "cheminformatics" 

## Create connection

## Check if connection already exists 

bq show --location=US --format=prettyjson --connection "cheminformatics-connection" > /dev/null 2>&1

status=$?

## if the connection exists, continue. otherwise, create the connection

if [ $status -eq 0 ]
 then
    echo "Connection cheminformatics-connect already exists"
 else
    echo "Creating connection cheminformatics-connect"
    bq mk --connection --display_name="Cheminformatics Connection" --connection_type=CLOUD_RESOURCE --location=US "cheminformatics-connection"
fi

## Get Service Account associated to Connection

SERVICE_ACCOUNT=$(bq show --location=US --format=prettyjson --connection "cheminformatics-connection" | jq -r '.cloudResource.serviceAccountId')

echo "Connection cheminformatics-connect service account: ${SERVICE_ACCOUNT}"

## Give service account the cloud run invoker role (necessary for cloud functions gen2)

PROJ=$(gcloud config list --format 'value(core.project)')

####### Begin creating cloud functions  

PERM="roles/cloudfunctions.invoker"

TIMEOUT=600s
MEMORY=512MB
MAX_INSTANCES=100

## start installation by folder

cd rdkit 

## install rdkit-fingerprint

gcloud beta functions deploy rdkit-fingerprint \
     --quiet --gen2 --region "us-east1" --entry-point rdkit_fingerprint --runtime python39 --trigger-http \
     --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES  \
     --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=rdkit

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-fingerprint --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-fingerprint" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

gcloud run services add-iam-policy-binding "rdkit-fingerprint" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics.rdkit_fingerprint(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 2500)'

## install rdkit-substructure-match

gcloud beta functions deploy rdkit-substructure-match \
    --quiet --gen2 --region "us-east1" --entry-point rdkit_substructure_match --runtime python39 --trigger-http \
    --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES  \
    --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=rdkit

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-substructure-match --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-substructure-match" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

gcloud run services add-iam-policy-binding "rdkit-substructure-match" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics.rdkit_substructure_match(fragment_smiles STRING, smiles STRING) RETURNS BOOL REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 2500)'

## install rdkit-molecular-descriptors

gcloud beta functions deploy rdkit-molecular-descriptors \
    --quiet --gen2 --region "us-east1" --entry-point rdkit_molecular_descriptors --runtime python39 --trigger-http \
    --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES  \
    --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=rdkit  

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-molecular-descriptors --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-molecular-descriptors" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

gcloud run services add-iam-policy-binding "rdkit-molecular-descriptors" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics.rdkit_molecular_descriptors(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 2500)'

## install rdkit-draw-svg

gcloud beta functions deploy rdkit-draw-svg \
     --quiet --gen2 --region "us-east1" --entry-point rdkit_draw_svg --runtime python39 --trigger-http \
     --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES \
     --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=rdkit

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-draw-svg --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-draw-svg" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

gcloud run services add-iam-policy-binding "rdkit-draw-svg" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics.rdkit_draw_svg(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 1000)'

## install rdkit-smiles-to-inchi

gcloud beta functions deploy rdkit-smiles-to-inchi \
     --quiet --gen2 --region "us-east1" --entry-point rdkit_smiles_to_inchi --runtime python39 --trigger-http \
     --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES \
     --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=rdkit

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-smiles-to-inchi --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-smiles-to-inchi" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

gcloud run services add-iam-policy-binding "rdkit-smiles-to-inchi" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics.rdkit_smiles_to_inchi(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 1000)'

## install rdkit-qed

gcloud beta functions deploy rdkit-qed \
    --quiet --gen2 --region "us-east1" --entry-point rdkit_qed --runtime python39 --trigger-http \
    --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES \
    --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=rdkit

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-qed --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-qed" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

gcloud run services add-iam-policy-binding "rdkit-qed" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics.rdkit_qed(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 1000)'


## move into biopython folder

cd ../biopython/

## install biopython-sequence-complement

gcloud beta functions deploy biopython-sequence-complement \
     --quiet --gen2 --region "us-east1" --entry-point biopython_sequence_complement --runtime python39 --trigger-http \
     --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES \
     --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=biopython

CLOUD_TRIGGER_URL=$(gcloud beta functions describe biopython-sequence-complement --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "biopython-sequence-complement" --gen2 --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

gcloud run services add-iam-policy-binding "biopython-sequence-complement" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.biopython_sequence_complement(sequence STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows=2500)'

## install biopython-sequence-translate

gcloud beta functions deploy biopython-sequence-translate \
      --quiet --gen2 --region "us-east1" --entry-point biopython_sequence_translate --runtime python39 --trigger-http \
      --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES \
      --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=biopython

CLOUD_TRIGGER_URL=$(gcloud beta functions describe biopython-sequence-translate --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "biopython-sequence-translate" --gen2 --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

gcloud run services add-iam-policy-binding "biopython-sequence-translate" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.biopython_sequence_translate(sequence STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows=2500)'

## install biopython-sequence-translate-to-stop

gcloud beta functions deploy biopython-sequence-translate-to-stop \
     --quiet --gen2 --region "us-east1" --entry-point biopython_sequence_translate_to_stop --runtime python39 --trigger-http \
     --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES \
     --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=biopython


CLOUD_TRIGGER_URL=$(gcloud beta functions describe biopython-sequence-translate-to-stop --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "biopython-sequence-translate-to-stop" --gen2 --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

gcloud run services add-iam-policy-binding "biopython-sequence-translate-to-stop" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.biopython_sequence_translate_to_stop(sequence STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows=2500)'

cd ..


## wait one minute for permissions to propagate
echo "Waiting for permissions to propagate ..."
sleep 90

