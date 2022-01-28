
## Summary

Combines hits from plasmidfinder (via [abricate](https://github.com/tseemann/abricate)), card database ARGs (via [abricate](https://github.com/tseemann/abricate)), [MOBsuite](https://github.com/phac-nml/mob-suite), and [Mobile Element Finder](https://pypi.org/project/MobileElementFinder/) to find MGEs with associated ARGs, and to isolate mobile plasmids with antibiotic resistance profiles. Provides diagrams of mobile plasmids with ARGs found in more than one input faster file. 


## Usage
Example run:

```bash hgt_pipeline.sh path/to/input/fasta/folder path/to/existing/output/directory 0.8 80 80 80 80 -5 0.00001 ```


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

Abricate runs automatically with e-value = 1e-20, and Mobile Element Finder runs automatically with %identity = 80%. 

## Output
The pipeline outputs several compiled results files as csv's, as well as preliminary visualizations of mobile plasmids with resistance found in more than one strain. 

1. all_hits.csv: Every MGE and ARG found


| Column      | Description |
| ----------- | ----------- |
| file_id      | Name of .fasta file       |
| contig_id   | contig in .fasta file where hit was found        |
| type   | ARG, type of MGE         |
| name   | ARG/MGE name        |
| src   | tool or database        |
| start   | start locus, if reported        |
| end   | end locus, if reported        |
| strand   | strand, if reported        |
| inc_rep   | plasmids only: inclusion group/replicon type , if reported        |
| mobility   | plasmids only: predicted mobility, if reported        |
| resistance   | ARGs only: drug class resistance conferred        |
| e_val   | e-value, if reported        |
| identity   | %identity, if reported        |
| coverage   | coverage, if reported        |

2. all_mges.csv: Every MGE found, with associated ARGs


| Column      | Description |
| ----------- | ----------- |
| file_id      | Name of .fasta file       |
| contig_id   | contig where hit was found        |
| type   | type of MGE         |
| name   | MGE name        |
| src   | tool or database        |
| start   | start locus, if reported        |
| end   | end locus, if reported        |
| strand   | strand, if reported        |
| inc_rep   | plasmids only: inclusion group/replicon type , if reported        |
| mobility   | plasmids only: predicted mobility, if reported        |
| ARGs   | all associated ARGs (plasmids: all ARGs on contig. Otherwise: all ARGs within 31kb)      |
| resistance   | all drug class resistance conferred by associated ARGs      |
| e_val   | e-value, if reported        |
| identity   | %identity, if reported        |
| coverage   | coverage, if reported        |

4. all_plasmids.csv: All plasmids found, with associated other MGEs, ARGs, and resistance]


| Column      | Description |
| ----------- | ----------- |
| file_id      | Name of .fasta file       |
| contig_id   | contig where plasmid was found        |
| plasmid   | name of plasmid         |
| inc_rep   |  inclusion group/replicon type , if reported        |
| ARGs   | all ARGs on contig      |
| resistance   | all resistance conferred by ARGs on contig    |
| mobility   |  predicted mobility, if reported        |
| MGEs   |  list of other MGEs on contig and location   |
| src   | plasmidfinder or MOB        |

5. plasmids_by_strain: One row per plasmid found in each strain, with contigs for that plasmid, and ARGs and resistance from all associated contigs in the strain


| Column      | Description |
| ----------- | ----------- |
| file_id      | Name of .fasta file       |
| plasmid   | name of plasmid         |
| inc_rep   | inclusion group/replicon type , if reported        |
| ARGs   | all ARGs on associated contigs      |
| resistance   | all resistance conferred by ARGs on associated contigs    |
| mobility   |  predicted mobility, if reported        |
| contig_id   | list of contigs in this strain identified as this plasmid   |
| src   | plasmidfinder or MOB        |


6. mobile_resistant_plasmids.csv: Summary file containing one row per plasmid name, with the number of strains found to contain that plasmid, the 5 most common drug class resistance conferred, and a list of the strains containing the plasmids. Only plasmids with associated ARGs and a predicted mobility != 'non-mobilizable' are included in the summary. 

| Column      | Description |
| ----------- | ----------- |
| plasmid   | name of plasmid   |
| num. strains   |   the number of strains containing the plasmid    |
| resistance   | 5 most common classes of drug resistance conferred by this plasmid |
| file_id   |   list of file names of strains containing the plasmid  |
