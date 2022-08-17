# cheminformatics-bq

Note: This is not an official Google product.

To install the cheminformatics package, open a cloud shell in the context of a project. 

Enable the following APIs

- gcloud services enable bigqueryconnection.googleapis.com
- gcloud services enable cloudfunctions.googleapis.com
- gcloud services enable run.googleapis.com
- gcloud services enable artifactregistry.googleapis.com
- gcloud services enable cloudbuild.googleapis.com

Installing the package requires the following permissions:

- BigQuery Admin
- Cloud Functions Admin
- Service Account User

As a user with above permissions, run **install_cheminformatics.sh**

This will:

- create a new dataset named "cheminformatics" 
- create a connection
- create the cloud functions (that wrap RDKit and Biopython)
- assign IAM permissions for the connection service account to invoke each cloud function
- assign IAM permissions for the connection service account to the underlying cloud run service (Cloud Functions Gen2)
- create the BQ functions  

Using the package requires the following permissions:

- BigQuery User
- BigQuery Connection User
- BigQuery Metadata Viewer
