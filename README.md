# cheminformatics-bq

To install the cheminformatics package, open a cloud shell in the context of a project. 

Then run **install_cheminformatics.sh**

This will:

- create a new dataset named "cheminformatics" 
- create a connection
- assign IAM permissions for the connection service account to cloud run (necessary for gen2 cloud functions)
- create the cloud functions (that wrap RDKit and Biopython)
- assign IAM permissions for the connection service account to invoke each cloud function
- create the BQ functions  

