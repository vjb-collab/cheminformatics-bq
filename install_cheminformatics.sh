#!/bin/bash 

## Create New Dataset to store the cheminformatics functions

bq mk --description "Dataset that will contain the cheminformatics functions" --dataset "cheminformatics" 

## Create Connection

bq mk --connection --display_name="Cheminformatics Connection" --connection_type=CLOUD_RESOURCE --location=US "cheminformatics-connection"

## Get Service Account associated to Connection

SERVICE_ACCOUNT=$(bq show --location=US --format=prettyjson --connection "cheminformatics-connection" | jq -r '.cloudResource.serviceAccountId')

echo "Connection created with service account: ${SERVICE_ACCOUNT}"

## Give service account the cloud run invoker role (necessary for cloud functions gen2)

PROJ=$(gcloud config list --format 'value(core.project)')

gcloud projects add-iam-policy-binding $PROJ --quiet --member=serviceAccount:$SERVICE_ACCOUNT --role=roles/run.invoker

## Create Cloud functions 

PERM="roles/cloudfunctions.invoker"

cd rdkit

## install rdkit-pattern-fingerprint

gcloud beta functions deploy rdkit-pattern-fingerprint --gen2 --region "us-east1" --entry-point rdkit_pattern_fingerprint --runtime python39 \
    --trigger-http --quiet --memory=512MB --timeout=240s --max-instances=9000  > /dev/null

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-pattern-fingerprint --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-pattern-fingerprint" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_pattern_fingerprint(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url)'

## install rdkit-morgan-fingerprint

gcloud beta functions deploy rdkit-morgan-fingerprint --gen2 --region "us-east1" --entry-point rdkit_morgan_fingerprint --runtime python39 \
     --trigger-http --quiet --memory=512MB --timeout=240s --max-instances=9000 > /dev/null

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-morgan-fingerprint --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-morgan-fingerprint" --gen2 --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_morgan_fingerprint(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url)'

# ## install rdkit-draw-svg

# gcloud functions deploy rdkit-draw-svg --entry-point rdkit_draw_svg --runtime python39 \
#     --trigger-http --quiet --memory=512MB --timeout=240s > /dev/null

# CLOUD_TRIGGER_URL=$(gcloud functions describe rdkit-draw-svg --format=json | jq -r '.httpsTrigger.url')

# gcloud functions add-iam-policy-binding "rdkit-draw-svg" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_draw_svg(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 50)'

# ## install rdkit-generate-maccs-keys

# gcloud functions deploy rdkit-generate-maccs-keys --entry-point rdkit_generate_MACCS_keys --runtime python39 \
#     --trigger-http --quiet --memory=512MB --timeout=240s > /dev/null

# CLOUD_TRIGGER_URL=$(gcloud functions describe rdkit-generate-maccs-keys --format=json | jq -r '.httpsTrigger.url')

# gcloud functions add-iam-policy-binding "rdkit-generate-maccs-keys" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_generate_maccs_keys(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url)'

# ## install rdkit-smiles-to-inchi

# gcloud functions deploy rdkit-smiles-to-inchi --entry-point rdkit_smiles_to_inchi --runtime python39 --trigger-http --quiet --memory=512MB --timeout=240s > /dev/null

# CLOUD_TRIGGER_URL=$(gcloud functions describe rdkit-smiles-to-inchi --format=json | jq -r '.httpsTrigger.url')

# gcloud functions add-iam-policy-binding "rdkit-smiles-to-inchi" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_smiles_to_inchi(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url)'

# ## install rdkit-descriptor-generic 

# gcloud functions deploy rdkit-descriptor-generic --entry-point rdkit_descriptor_generic --runtime python39 \
#     --trigger-http --quiet --memory=512MB --timeout=240s

# CLOUD_TRIGGER_URL=$(gcloud functions describe rdkit-descriptor-generic --format=json | jq -r '.httpsTrigger.url')

# gcloud functions add-iam-policy-binding "rdkit-descriptor-generic" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

# ## install descriptor library

# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_descriptor_exact_mol_wt(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, user_defined_context = [("rdkit-function", "ExactMolWt")] )'
# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_descriptor_mol_log_p(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, user_defined_context = [("rdkit-function", "MolLogP")] )'
# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_descriptor_fraction_csp3(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, user_defined_context = [("rdkit-function", "FractionCSP3")] )'
# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_descriptor_num_aliphatic_rings(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, user_defined_context = [("rdkit-function", "NumAliphaticRings")] )'
# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_descriptor_heavy_atom_count(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, user_defined_context = [("rdkit-function", "HeavyAtomCount")] )'

# ## install rdkit qed

# gcloud functions deploy rdkit-qed --entry-point rdkit_qed --runtime python39 --trigger-http --quiet --memory=512MB --timeout=240s > /dev/null

# CLOUD_TRIGGER_URL=$(gcloud functions describe rdkit-qed --format=json | jq -r '.httpsTrigger.url')

# gcloud functions add-iam-policy-binding "rdkit-qed" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.rdkit_qed(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url)'

# ## move into biopython folder

# cd ../biopython/

# ## install biopython-sequence-complement

# gcloud functions deploy biopython-sequence-complement --entry-point biopython_sequence_complement --runtime python39 --trigger-http --quiet --memory=512MB --timeout=240s > /dev/null

# CLOUD_TRIGGER_URL=$(gcloud functions describe biopython-sequence-complement --format=json | jq -r '.httpsTrigger.url')

# gcloud functions add-iam-policy-binding "biopython-sequence-complement" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} 

# bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE FUNCTION cheminformatics.biopython_sequence_complement(sequence STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url)'

# ## leave biopython 

# cd .. 

## wait one minute for permissions to propagate
echo "Waiting for permissions to propagate ..."
sleep 90

# 100 million in 16 minutes
# 978 million in 1 hour

