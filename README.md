#Summary
Combines hits from plasmidfinder (via abricate), card database ARGs (via abricate), mob_recon, and Mobile Element Finder to find MGEs with associated ARGs, and to isolate mobile plasmids with antibiotic resistance profiles. Provides diagrams of mobile plasmids with ARGs found in more than one input faster file. 

#Usage
Example run:
```bash hgt_pipeline.sh path/to/input/fasta/folder path/to/existing/output/directory 0.8 80 80 80 80 -5 0.00001 1```


Command line args:
1. Location of fasta files
2. results directory
3. mefinder coverage: Float, 0 to 1
4. Abricate coverage: Int, 0 to 100
5. Mob coverage: Int, 0 to 100
6. Abricate percent identity: Int, 0 to 100
7. Mob percent identity: Int, 0 to 100
8. Mefinder e-value: Int, 1eX, example -6 is 1e-6
9. Mob evalue: Float, min allowed 0.00001
10. Visuals: 0 = no, 1 = yes

Abricate runs automatically with e-value = 1e-20, and Mobile Element Finder runs automatically with %identity = 80%. 
