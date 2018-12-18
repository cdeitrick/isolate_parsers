# isolate_parsers

## Usage
```
python isolateset_parser.py [-h] [-i FOLDER] [--no-fasta] [-w WHITELIST]
                            [-b BLACKLIST] [-m SAMPLE_MAP] [--filter-1000bp]

optional arguments:
  -h, --help            show this help message and exit
  -i FOLDER, --input FOLDER
                        The breseq folder to parse.
  --no-fasta            Whether to generate an aligned fasta file of all snps
                        in the breseq VCF file.
  -w WHITELIST, --whitelist WHITELIST
                        Samples not in the whitelist are ignored. Either a
                        comma-separated list of sample ids for a file with
                        each sample id occupying a single line.
  -b BLACKLIST, --blacklist BLACKLIST
                        Samples to ignore. See `--whitelist` for possible
                        input formats.
  -m SAMPLE_MAP, --sample-map SAMPLE_MAP
                        A file mapping sample ids to sample names. Use if the
                        subfolders in the breseqset folder are named
                        differently from the sample names. The file should
                        have two columns: `sampleId` and `sampleName`,
                        separated by a tab character.
  --filter-1000bp       Whether to filter out variants that occur within
                        1000bp of each other. Usually indicates a mapping
                        error.
```

## Input
The scripts expect a folder of individual breseq runs, with each folder named after the isolate/sample.
The scipts only require the `output.vcf`, `annotated.gd`, and `index.html` files located in each folder.
Example folder:
```
    .breseq_folder
    |-- sample1
    |   |-- data
    |   |   |-- output.vcf
    |   |-- output
    |   |   |-- index.html
    |   |   |-- evidence
    |   |   |   |-- annotated.gd
    |-- sample2
    |   |-- data
    |   |   |-- output.vcf
    |   |-- output
    |   |   |-- index.html
    |   |   |-- evidence
    |   |   |   |-- annotated.gd
    |-- sample3
    |   |-- data
    |   |   |-- output.vcf
    |   |-- output
    |   |   |-- index.html
    |   |   |-- evidence
    |   |   |   |-- annotated.gd
```

## Output
The scripts generate an excel file in the breseq run folder with 4 sheets: `comparison`, `variant`, `coverage`, and `junction`.
The `variant`, `coverage`, and `junction` tables are just the concatenated tables of all samples in the breseq run.

### Comparision table
A table in which every row represents a single mutation seen in the sample callset 
and samples are represented by columns with the alternate sequence for each sample.

| Sample1 | Sample2 | Sample3 | annotation                | description                                                                                     | gene                    | locusTag          | mutationCategory  | position | presentIn | presentInAllSamples | ref | seq id    | 
|---------|---------|---------|---------------------------|-------------------------------------------------------------------------------------------------|-------------------------|-------------------|-------------------|----------|-----------|---------------------|-----|-----------| 
| GG      | GG      | GG      | intergenic (+65/+20)      | putative lipoprotein/putative hydrolase                                                         | PFLU0045 - / - PFLU0046 | PFLU0045/PFLU0046 | small_indel       | 45881    | 3         | 1                   | G   | NC_012660 | 
| CC      | CC      | CC      | intergenic (+17/-136)     | microcin-processing peptidase 1. Unknown type peptidase. MEROPS family U62/hypothetical protein | PFLU0872 - / - PFLU0873 | PFLU0872/PFLU0873 | small_indel       | 985333   | 3         | 1                   | C   | NC_012660 | 
|         |         |         | intergenic (+57/+21)      | hypothetical protein/putative helicase                                                          | PFLU3154 - / - PFLU3155 | PFLU3154/PFLU3155 | small_indel       | 3447986  | 3         | 1                   |     | NC_012660 | 
| A       | A       | G       | M350I (ATG-ATA)           | putative GGDEF domain signaling protein                                                         | PFLU3571 -              | PFLU3571          | snp_nonsynonymous | 3959631  | 2         | 0                   | G   | NC_012660 | 
| A       | A       | C       | T238P (ACC-CCC)           | hybrid sensory histidine kinase in two-component regulatory system with UvrY                    | PFLU3777 -              | PFLU3777          | snp_nonsynonymous | 4173231  | 1         | 0                   | A   | NC_012660 | 
| G       | G       | GG      | coding (322/1476 nt)      | putative two-component system response regulator nitrogen regulation protein NR(I)              | PFLU4443 -              | PFLU4443          | small_indel       | 4908233  | 1         | 0                   | G   | NC_012660 | 

### Aligned fasta files

The scripts also generates 3 fasta files (`breseq.snp.fasta`, `breseq.amino.fasta`, `breseq.codon.fasta`)
with all nonsynonymous snps from each sample represented by the replacement bases, amino acids, and codons.
Example:
```
>reference
GA
>Sample1
AA
>Sample2
AA
>Sample3
GC
```