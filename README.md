# ======================================================================
# iStrength.R (v11) README
# ======================================================================

iStrength.R is an R script to compute chromatin interaction enrichment metrics from Hi-C data at genomic domains (TADs, loops, anchors).

It quantifies:
- Boundary strength (iStrength)
- Intra-domain interactions
- Inter-domain interactions
- Stripe patterns (upstream/downstream)
- Loop intensities
- Multiple normalization strategies (obs, OE, log, median)

# ======================================================================
# 1. QUICK OVERVIEW (WHAT IT ANALYZES)
# ======================================================================

Genome view:

        Flanking upstream                 TAD / DOMAIN                  Flanking downstream
<---- window ---->      |<------ domain ------>|      <---- window ---->

Upstream stripe        Boundary START      INTRA-DOMAIN       Boundary END       Downstream stripe
====****_____          |================|================|          _____****====

                       INTER-DOMAIN (background / insulation)
---------------------- ---------------- ---------------- ----------------------

Loop view:

[ Anchor A ] ======== loop interaction ======== [ Anchor B ]

Legend:
*  = stripe enrichment signal
=  = intra-domain interaction
-  = inter-domain/background signal
Anchors = loop endpoints (inside or boundaries)

# ======================================================================
# 2. HOW TO RUN
# ======================================================================

Rscript iStrength.v11.R config.yaml experiment_name

Example:
Rscript iStrength.v11.R config.yaml exp1

Where:
- config.yaml = configuration file
- experiment_name = key inside config under "experiments"

# ======================================================================
# 3. EXAMPLE CONFIGURATION (config.yaml)
# ======================================================================

general:
  cores: 8
  genome: mm10
  hic_directory: /path/hic/
  regions_directory: /path/regions/
  resolution: 10000
  pairfile: bedpe
  window: 200000
  window_ratio: 2
  dmin: 50000
  stripe_width: 50000
  skiping: 0
  regions_header: TRUE
  keep_id: TRUE

experiments:
  exp1:
    hic_file: sample.hic
    regions: loops.bedpe
    output_name: test_run

# ======================================================================
# 4. INPUTS
# ======================================================================

Hi-C file (.hic)
- Extracted using straw:
  - observed counts
  - observed/expected (OE)
  - raw counts

Regions file:
Defines anchors (BED or BEDPE)

BED:
chr   start   end   id

BEDPE:
chr1 start1 end1   chr2 start2 end2   id

# ======================================================================
# 5. PIPELINE (STEP BY STEP)
# ======================================================================

1. LOAD DATA
   - config.yaml
   - genome sizes
   - Hi-C matrices
   - genomic regions

2. FILTER REGIONS
   - remove small domains (< dmin)
   - define analysis window

3. EXTRACT Hi-C SIGNAL
   - observed
   - OE
   - raw counts

4. COMPUTE FEATURES

   TAD interior:
   =====================

   Inter-domain signal:
   ---------------------

   Boundary strength:
   iStrength = (intra - inter) / (intra + inter)

   Stripe signal:
   ====****____ (upstream)
   ____****==== (downstream)

   Loop signal:
   [A] ======= [B]

5. BACKGROUND NORMALIZATION
   - remove stripes
   - remove loops
   - compute baseline domain signal

6. AGGREGATE RESULTS
   - per chromosome
   - per domain

7. EXPORT OUTPUTS

# ======================================================================
# 6. OUTPUT FILES
# ======================================================================

<iS>.bed
- full domain metrics

<iS>.bedpe
- paired anchor representation

<iS>_summary.bed
- simplified boundary-focused metrics

# ======================================================================
# 7. KEY METRICS
# ======================================================================

Boundary strength:
- iStrength_start_obs
- iStrength_end_obs
- iStrength_*_oe
- iStrength_*_log

Loop metrics:
- loop_obs
- loop_oe
- loop_ratio

Stripe metrics:
- up_stripe
- down_stripe
- stripe_ratio

TAD metrics:
- tad_obs
- tad_log
- tad_oe

# ======================================================================
# 8. PARALLELIZATION
# ======================================================================

- One chromosome per core (mclapply)
- Higher resolution = slower runtime
- straw extraction is main bottleneck

# ======================================================================
# AUTHOR
# ======================================================================

Oscar Amaury Aguilar Lomas
2025
