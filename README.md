# ROBOINDEX

General developmental documentation for the RoboIndex pipeline.

---

## Running the pipeline

Refer to this page for general guidance on how to run the pipeline: [How to run ROBOINDEX](https://stowersinstitute.sharepoint.com/sites/ComputationalBiologyPublic/SitePages/RoboIndex.aspx)

---

## Code Structure

The source code for PRIME is stored under: /n/ngs/tools/ROBOINDEX/

In this and most nextflow pipelines:
- main.nf defines the general logic
- workflows/ contain all workflows defined in main.nf
- modules/ contain all processes defined in main and workflows under workflows/
- bin/ contain all the python and r scripts used in processes under module/ and workflows/
- nextflow.config contains all process parameters for slurm resource allocations 

---

## General Logic 

The parameters and the general logic of the pipeline is defined in "main.nf". 
In this file, I specified in a series of workflows to the type of index one may wish to build.

***Add additional workflow to include new order types further down the line.***

---

## Workflows

Workflows of different index types are designed to be independent of each other.
I would encourage people to add new workflows when adding new index types. 
Right now all the workflows are stored in Genome_Index.nf under the workflow directory. 

---

## Modules

Processes are grouped based on function. 
- core: Processes that are essential in Stower's application
- bacteria: Process used in building bacteria reference
- fetch-ensembl: processes used in downloading genomes and building references from ensembl 

---

## Run Environment

PRIME is running within a conda environment stored here: /home/compbio_svc/miniconda3/envs/RoboIndex_ENV

---
