# telomere_finder

Telomere finder in edges of contigs.

```
python3 telomere_finder.py -i contigs.fasta
```

It search **CCCTAAA** and its complement repeat sequences for plant.



optional arguments:
  -h, --help            show this help message and exit
  -i str                input fasta file of genome
  -e int                Length of the end that is skipped when non-telomere
                        sequences (default=100)
  -a int                Distance to allow not telomere sequences (default=28)
  -l int                Minimum contig length (default=100000)
  -t comma_separated_str
                        telomere sequence. comma separated.
                        (default=CCCTAAA,TTTAGGG for plant)
  -o str                output file name of bed format.

OUTPUT fasta_ID telomere_sequence F/R Dist_from_edge num_telomere


