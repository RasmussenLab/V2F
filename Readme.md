# Multimodal deep learning connects polygenic metabolic disease risk to adipocyte regulatory circuits

Code for the analyses in *Multimodal deep learning connects polygenic metabolic disease risk to
adipocyte regulatory circuits*, a variant-to-function (V2F) study of adipogenesis and metabolic
disease. Most scripts run in protected environments against data that is not openly available,
so this repository is primarily a durable record of the code rather than a runnable pipeline.

**The study.** Adipose-derived mesenchymal stem cells (AMSCs) from 37 donors, differentiating
towards mature adipocytes, profiled across 269 samples for gene expression, LipocyteProfiler
morphology, polygenic risk scores and cardiometabolic variants. The four analysis stages follow
one thread from that dataset to a mechanism:

1. **Integrate and generate hypotheses.** A variational autoencoder integrates the modalities;
   *in silico* perturbations link insulin-resistance polygenic risk to *SHC2* expression, and
   the *C5orf67* variant rs3843467 to *SHC2*, *MAP3K1* and *SETD9*.
2. **Prioritise variants.** Colocalization across four glycaemic traits and fine mapping of T2D
   and fasting insulin over *C5orf67* narrow the signal to rs256903 and a set of variants in
   near-complete linkage with it.
3. **Predict what the variants do.** Sequence-to-function models place them in open chromatin
   next to AP-1 and SP/KLF motifs, and map their effect on *MAP3K1* to adipose contexts.
4. **Validate independently.** A single-cell village of differentiating AMSCs confirms the motif
   enrichment and finds context-dependent eQTLs under glucose stress.

## Where to find each part of the paper

| Paper | Code |
| --- | --- |
| Fig. 1, Suppl. Figs. 1–3 — VAE integration and SHAP feature importance | [1_data_integration/](scripts/1_data_integration/) |
| Fig. 2a–g, Suppl. Figs. 4–5 — *SHC2* as a convergence point, genotype–expression follow-up | [1_data_integration/](scripts/1_data_integration/) |
| Fig. 2h–j, Suppl. Figs. 6–7 — colocalization and fine mapping | [2_colocalization_and_finemapping/](scripts/2_colocalization_and_finemapping/) |
| Fig. 3a — ChromBPNet accessibility and contribution scores | [3_.../ChromBPNet/](scripts/3_sequence_to_function_modeling/ChromBPNet/) |
| Fig. 3b, Fig. 4, Suppl. Figs. 8–17 and 21 — Borzoi and AlphaGenome, and their cross-comparison | [3_.../AlphaGenome_Borzoi/](scripts/3_sequence_to_function_modeling/AlphaGenome_Borzoi/) |
| Suppl. Figs. 18–20 — single-cell eQTLs and caQTLs | [4_adipogenic_village_analyses/](scripts/4_adipogenic_village_analyses/) |

## Repository layout

**[1_data_integration/](scripts/1_data_integration/)** —
[AMSC_MOVE.ipynb](scripts/1_data_integration/AMSC_MOVE.ipynb), the main notebook of the project:
preprocessing, integration with [MOVE](https://www.nature.com/articles/s41587-022-01520-x), SHAP
feature importance, the *in silico* perturbation analyses, and the downstream *SHC2* / *MAP3K1*
genotype–expression analyses.

**[2_colocalization_and_finemapping/](scripts/2_colocalization_and_finemapping/)** — each stage
reads from `dat/` and writes to its own `interim/` and `results/`, none of which are tracked.

- **2a_colocalization** —
  [SharePro_signed.ipynb](scripts/2_colocalization_and_finemapping/2a_colocalization/SharePro_signed.ipynb)
  builds the signed LD matrix, harmonises the four GWAS onto it, and runs
  [SharePro](https://doi.org/10.1093/bioinformatics/btae295) across T2D, fasting insulin, fasting
  glucose and 2-hour glucose (`res_signed` is the three-trait model, `res_signed_2h` adds 2-hour
  glucose).
- **2b_finemapping** —
  [CARMA_T2D_signed.R](scripts/2_colocalization_and_finemapping/2b_finemapping/CARMA_T2D_signed.R)
  fine-maps T2D and fasting insulin with
  [CARMA](https://www.nature.com/articles/s41588-023-01392-0);
  [CARMA_seed_stability.R](scripts/2_colocalization_and_finemapping/2b_finemapping/CARMA_seed_stability.R)
  reruns it across 20 seeds and reports which conclusions survive, since CARMA's stochastic
  search sets no seed of its own and this locus holds many perfectly collinear variants.
- **dat** — what 2a writes and 2b reads: per-trait effect sizes and standard errors, the signed
  European LD matrix around rs3843467, and `variant_order.tsv`, its row index, which records the
  allele every sign and every beta refers to. The matrix carries no labels of its own, so those
  two belong together. `variant_list.txt` is the only input fixed in advance. The matrix comes
  from `plink --r square` on 1000 Genomes EUR genotypes rather than from LDlink, whose endpoints
  return only unsigned R² or D′ — SharePro and CARMA need the *signed* correlation to tell
  whether two correlated effects reinforce or cancel.

**[3_sequence_to_function_modeling/](scripts/3_sequence_to_function_modeling/)**

- **ChromBPNet** —
  [ChromBPNet-C5orf67.ipynb](scripts/3_sequence_to_function_modeling/ChromBPNet/ChromBPNet-C5orf67.ipynb):
  bias-corrected ATAC-seq profile prediction on bulk visceral AMSCs, DeepLIFT/DeepSHAP
  contribution scores and *in silico* substitutions. The variant-centred regions and SNP list it
  takes as inputs are built in `auxilliary_files/`.
- **AlphaGenome_Borzoi** —
  [C5orf67_MAP3K1_Borzoi_and_AlphaGenome.ipynb](scripts/3_sequence_to_function_modeling/AlphaGenome_Borzoi/C5orf67_MAP3K1_Borzoi_and_AlphaGenome.ipynb),
  both long-range models in one kernel so their results can be compared directly: Borzoi
  predictions of *MAP3K1* across tracks with the hypergeometric test for adipose enrichment,
  AlphaGenome variant scoring, ISM and haplotype editing, and a cross-model section on why the
  two disagree on direction. The variant list is not typed in — it is derived at run time from
  the CARMA and SharePro outputs.

**[4_adipogenic_village_analyses/](scripts/4_adipogenic_village_analyses/)** —
[Adipogenic_village_analyses.v1.ipynb](scripts/4_adipogenic_village_analyses/Adipogenic_village_analyses.v1.ipynb):
the AMSC cell village 10X-Multiome analyses as presented in the manuscript — motif enrichment in
peaks upregulated in the adipogenic subpopulation, and the eQTL / caQTL analyses around
*C5orf67*.

**[configs/MOVE/](configs/MOVE/)** — the `data`, `experiment` and `tasks` YAML files used to run
the MOVE models.

## Running the analyses

This repository tracks the analysis code, the MOVE configuration and this Readme, and nothing
else. Data and outputs are git-ignored — the `dat/`, `data/` and `auxilliary_files/` folders each
stage reads from, and the `interim/` and `results/` folders it writes — so a fresh clone holds no
inputs, figures or tables. Each stage rebuilds what the one before it produced, and stage 3 reads
stage 2's output, so run them in order, or point a stage at a folder that already holds what it
needs.

Also kept out of version control, and so present only in a working copy: the manuscript, local
copies of the [SharePro](https://doi.org/10.1093/bioinformatics/btae295) and
[CARMA](https://www.nature.com/articles/s41588-023-01392-0) packages, the ChromBPNet tracks
exported for the [WashU Epigenome Browser](https://epigenomegateway.wustl.edu/) in Fig. 3a, and
the superseded results from before the LD matrix carried a sign.

The AlphaGenome notebook needs an API key. Put a line reading `ALPHA_GENOME_API_KEY=your-key` in a
`.env` file at the repository root; VS Code loads it into the kernel automatically, and `.env` is
git-ignored. The notebook's *Running locally* section covers the rest.

## Donor privacy and reproducibility

**No donor or patient identifiers belong in this repository.** The AMSC data is individual-level
and stays in the protected environment where the notebooks are run. The code refers to donors
only through placeholders such as `<DONOR_ID>` and through column names such as `Patient_ID`;
cloud paths are likewise reduced to placeholders. Everything tracked here is variant-level. This
matters most in `1_data_integration` and `4_adipogenic_village_analyses`, the two stages that
touch donor-level data directly.

`.gitignore` blocks identifier files, MOVE working directories, access-controlled source data,
API keys and all regenerable output. Before committing, check that nothing slipped through:

```bash
git add -An | sed 's/^add //' | tr -d "'" | grep -iE 'donor|patient|private|sumstat|\.vcf'
```

**Strip notebook outputs before committing** — they can carry identifiers, absolute paths and
intermediate data, and they make diffs unreadable:

```bash
jupyter nbconvert --clear-output --inplace path/to/notebook.ipynb
# or once, to have it happen on every commit:
pip install nbstripout && nbstripout --install
```

## Data availability

The AMSC data is not openly available. Transcriptomics from the LipocyteProfiler study are on GEO
under [GSE184089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184089), and the
high-content imaging data in the Cell Painting Gallery on the AWS Registry of Open Data under
`cpg0011`. GWAS summary statistics used for colocalization and fine mapping:

- Fasting insulin — [GCST90002238](https://www.ebi.ac.uk/gwas/studies/GCST90002238)
- Two-hour glucose — [GCST90002227](https://www.ebi.ac.uk/gwas/studies/GCST90002227)
- Fasting glucose — [GCST90002232](https://www.ebi.ac.uk/gwas/studies/GCST90002232)
- Type 2 diabetes — [DIAMANTE](https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_DIAMANTE_eu)

## Related code

MOVE is available as the `move-dl` pip package. Code and configuration for creating and running
MOVE models are also at
[RasmussenLab/VAEs_for_biomedical_data_integration](https://github.com/RasmussenLab/VAEs_for_biomedical_data_integration),
with a
[tutorial notebook](https://colab.research.google.com/github/RasmussenLab/VAEs_for_biomedical_data_integration/blob/main/scripts/Tutorial_VAEs_for_biomedical_data_integration.ipynb)
covering the main findings.
