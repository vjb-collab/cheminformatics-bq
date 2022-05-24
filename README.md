# cheminformatics-bq

To install the cheminformatics package, open a cloud shell in the context of a project. 

Then run **install_cheminformatics.sh**

This will:

- create a new dataset named "cheminformatics" 
- create a connection
- create the cloud functions (that wrap RDKit and Biopython)
- assign IAM permissions for the connection service account to invoke each cloud function
- create the BQ functions  

