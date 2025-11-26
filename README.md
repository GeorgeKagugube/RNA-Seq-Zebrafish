# Trascriptomic Profile of Mn overloaded Zerbafish Brain
End-to-end bulk RNA-seq pipeline from experimental design, RNA sample collection, and analysis pipeline, including raw FASTQ files QC to differential gene expression, GO term, and Pathway analysis.

<!-- Optional badges -->
<!--
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.3.0-blue)]()
[![Build Status](https://github.com/USER/REPO/actions/workflows/ci.yml/badge.svg)]()
-->

---

## ✨ Features

- ✅ First key feature (e.g. “End-to-end QC → alignment → DE → plots”)
- ✅ Second key feature
- ✅ Third key feature

---

## 📁 Project structure

Brief description of the main folders / files.

```text
PROJECT_NAME/
├─ data/              # Input data / placeholders / examples (no raw data if sensitive)
├─ results/           # Output tables, figures, logs
├─ src/               # Core code (R / Python / etc.)
├─ notebooks/         # Exploratory analyses / reports (Rmd, qmd, ipynb)
├─ config/            # Config files (YAML/JSON) for parameters, paths, etc.
├─ env/               # Conda environment or requirements files
└─ README.md          # This file


## Background
<div style='text-align: right;'>
<p style='text-align: right;'>The trace element Manganese (Mn) is involved in key biochemical reactions central to normal cellular and systemic function. 
Like other trace metals, its cellular and systemic concentration is tightly regulated, mainly through dietary absorption 
and excretion via GI Track and hepatobiliary routes. Mn-driven neurotoxicity is implicated in neurodegenerative diseases 
such as Parkinson's, but the mechanisms remain elusive. Mitochondrial dysfunction, disturbed calcium physiology, and 
oxidative stress are potential pathways involved in Mn overload. </p>

<p style='text-align: right;'>
The influx transporters SLC39A14/8 and efflux SLC30A10 have higher Mn transportation preferences and are central to Mn 
homeostasis. Using a loss of function mutation of slc39a14<sup>-/-</sup>
compared to its slc39a14<sup>+/+</sup> controls in a zebrafish, Tuschl et al. (2020) recreated Mn overload in the brain 
reported in clinical settings with a compromised movement disorder typical of Mn overload patients. Studying the response 
and sensitivity to the Mn challenge of these mutants creates a representative model organism for studying Mn-related 
neurotoxic mechanisms.  A zebrafish also presents a good model organism given its well-annotated genome, ease of imaging (allowing for a non-invasive look into its brain at single-cell resolution), and high reproductive and growth rates.
</p>

**Aim:** To identify the key molecular players and pathways in Mn overload-driven neurotoxicity by leveraging transcriptomics and computational methods. 

The pipelines shared here were used to analyse the transcriptome landscape of slc39a14<sup>-/-</sup>and slc39a14<sup>+/+</sup> in exposed compared to unexposed conditions. 

</div>

## Project Pipeline
# ![Workflow](https://github.com/GeorgeKagugube/bulk-rnaseq-zebrafish-mouse-human/blob/main/images/RNA%20ANALYSIS%20WORKFLOW.jpeg)
