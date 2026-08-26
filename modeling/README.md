## Folder for predictions from state-of-the-art models

### Enformer and AlphaGenome
- Enformer and Alphagenome predictions were taken as described by the original publication and the corresponding github repositories.
- The biggest models were used and the mean and max absolute difference for all the tracks was used as variant score
- Example script is planed to be published soon

### CADD v1.6 and v1.7
- For CADD the webserver was used to obtain variant effect predictions and the PHRED scaled scores were used


### Fine-tuning AlphaGenome
- For fine-tuning AlphaGenome the blog post and the respective github repository was used and the encoder was fine-tuned in a 3-fold cross validation manner on the MPRA data using the minimal promoter sequence and augmentation techniques similar to those described in the blog post. Example code requires own AlphaGenome API key and is planed to be published soon.
