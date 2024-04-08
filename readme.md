# bioTools

## What is it

A minimal C++ lib for reading / writing / editing common bioinformatics formats, supporting :

- bed
- vcf
- fasta
- wiggle

Implements also a new one, AOE for area of effect, a 5-col bed + a 6th one for the zero of the interval

Also builds two small softs :

- countFeatures counts occurences of bed entries along aoe
- intersectKeepingNames fills a void in bedtools intersect, allowing to intersect two beds, keeping the A entry as it is except for its start and end

## How to build

```
cmake .
cmake --build .
```