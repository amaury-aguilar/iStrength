# iStrength.R

## Quick Description

iStrength.R is an R script that computes chromatin interaction strength across genomic domains using Hi-C data. It quantifies how strongly regions such as TAD-like domains, loop anchors, and stripe structures interact internally compared to their surrounding genomic background, producing an interaction strength metric called iStrength that reflects domain insulation and boundary strength.

## How to Use

Run:
Rscript iStrength.R config.yaml experiment_name

Example:
Rscript iStrength.R config.yaml WT_sample

## Input

A YAML configuration file defines parameters and input paths.

general:
  cores: 8
  genome: mm10
  hic_directory: /path/to/hic/
  regions_directory: /path/to/regions/
  resolution: 10000
  pairfile: bedpe
  choose_window: fixed
  window: 200000
  window_ratio: 2
  dmin: 50000
  window_anchors_in: 50000
  window_anchors_out: 50000
  stripe_width: 50000
  skiping: 0
  regions_header: TRUE
  keep_id: TRUE

experiments:
  WT_sample:
    hic_file: sample.hic
    regions: regions.bedpe
    output_name: WT_sample

## Output

Results are written to:
iStrength/<experiment_name>/<resolution>/

Outputs include:
- .bed (full per-region statistics)
- .bedpe (domain representation in paired format)
- .iStrength (boundary strength summary)

## What the script does

The script evaluates chromatin organization using Hi-C contact maps by measuring intra-domain, inter-domain, boundary, stripe, and loop interaction signals.

UPSTREAM        DOMAIN              DOWNSTREAM
====*****====|===========|====*****====
   signal         core        signal

For each genomic region it computes:
- intra-domain interaction signal
- inter-domain interaction signal
- boundary enrichment (start and end)
- stripe enrichment
- loop interaction enrichment

## Boundary strength (iStrength)

START BOUNDARY:
====*****====|

END BOUNDARY:
|====*****====

iStrength = (intra - inter) / (intra + inter)

This measures how strongly a domain is insulated from its surrounding chromatin environment.

## Stripes

====*****----------
      \\\\\\\\\\\\
       \\\\\\\\\\\\

Directional interaction enrichment near domain boundaries.

## Loops

ANCHOR A ================= ANCHOR B
          \\\\\\\\\\\\\\\\\

Loop enrichment is computed using interaction windows around anchor pairs.

## Workflow

CONFIG FILE → LOAD REGIONS → LOAD Hi-C MATRIX (strawr) → ITERATE CHROMOSOMES → ITERATE REGIONS → COMPUTE METRICS → MERGE RESULTS → WRITE OUTPUT FILES

## Key concept

iStrength captures the balance between internal and external chromatin interactions, providing a quantitative measure of domain insulation and boundary strength in Hi-C data.
