# DESeq for Metagenome
The Taxamerg.biom file was created by merging the biom files from the individual samples, which were downloaded from MG-RAST. To merge the different samples biom files use  [MG-RAST tools](https://github.com/MG-RAST/MG-RAST-Tools)  (file 'mg-biom-merge.py' in 'scripts').

The RottnestFebMap.txt file is the map (metadata) for the samples data set. It has to contain the samples ID, which are given by MG_RAST (fist collumn in the file), a code for the samples, which is the name I gave to my samples (second collumn in the file) and the other collumns are the factors that can be tested.
