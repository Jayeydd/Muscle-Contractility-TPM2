# Muscle-Contractility-TPM2
## Molecular Mechanisms of Muscle Contractility: Protein Structure and Mutation Analysis
**Course:** Cell and Molecular Biology

**Activity:** Individual Laboratory

**Author:** Jade Angela Suan

**Date Accessed (all databases):** 2026-10-28

---

## (a) Assignment
| Item | Detail |
|---|---|
| Assignment number | 18 |
| Gene | TPM2 |
| NCBI Gene ID | 7169 (chromosome 9, cytogenetic band 9p13.3) |
| Transcript | NM_003289.4 |
| Matching protein | NP_003280.2 (tropomyosin beta chain, isoform Tpm2.2st), 284 aa |
| Variant | NM_003289.4(TPM2):c.397C>T (p.Arg133Trp), missense |
| ClinVar Variation ID | 12463 |
| UniProt | P07951 |

---

## (b) Research Question
How does the missense variant p.Arg133Trp in beta-tropomyosin relate to the regulation of skeletal muscle contraction, and what does the evidence show versus what is only predicted?

---

## (c) Methods
Software: R / RStudio. 
Packages: `rentrez`, `ggplot2`. 
All code is in `analysis.R`.

1. **Project setup** — Created folder `Muscle_Contractility_TPM2` with `data/`, `results/`, `figures/` and `analysis.R`.
```r
# Run installation only once on your computer:
install.packages(c("rentrez", "ggplot2"))
library(rentrez)
library(ggplot2)
# CHANGE THESE FOUR VALUES to your assignment (example: ACTA1):
gene_symbol <- "TPM2"
transcript_accession <- "NM_003289.4"
clinvar_id <- "12463"
variant_protein <- "p.Arg133Trp"
# Use "" for splice variants for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}

```
2. **Gene record**. Searched NCBI Gene with `entrez_search(db = "gene", term = "TPM2[Gene Name] AND txid9606[Organism]")`, checked the exact symbol with `entrez_summary()`, and saved `results/gene_summary.txt.`

3. **Protein from the assigned transcript.** Fetched NM_003289.4 with `entrez_fetch(db = "nuccore", rettype = "gb")`, extracted the `/protein_id` (NP_003280.2) with `regmatches()`, then fetched the FASTA with `entrez_fetch(db = "protein", rettype = "fasta")`. Saved `data/TPM2_protein.fasta` and `results/protein_accession.txt`.

4. **Composition.** Read the FASTA, counted the 20 standard amino acids with `table(factor(aa, levels = standard))`, saved `results/amino_acid_composition.csv` and plotted with `ggplot2::geom_col()` to `figures/01_amino_acid_composition.png`.

5. **ClinVar review.** Opened https://www.ncbi.nlm.nih.gov/clinvar/variation/12463/ and recorded the variant name, HGVS, type, conditions, classification, review status, submissions and citations (section e).

6. **Variant position figure.** Extracted the residue number from `variant_protein` with `regmatches()` and plotted it on the 284-aa line with `ggplot2` (`geom_segment, geom_point, annotate`), saved with `ggsave(..., bg = "white")` as `figures/02_variant_position.png`.

7. **Structure figure.** Opened the AlphaFold DB model AF-P07951-F1 (UniProt P07951), highlighted residue 133 (viewer popup read "P07951: 133"), and saved a screenshot as `figures/protein_structure_or_splice_diagram.png`.

**Errors corrected:** The PDF’s code lost its line breaks when copied, causing errors. Three variable assignments (`transcript_accession, clinvar_id, variant_protein`) were on one line, so R could not read `transcript_accession` when Step 3 tried to fetch the transcript. I fixed this by re-declaring the four variables on separate lines before Step 3 (and in the final script, at the top). Other problems were `gene_id <- NA_character_` inside a comment, the `if (nzchar(variant_protein)) {` line commented out, and an extra quote in `ggsave`. I put each statement on its own line and fixed the quote. I also added `bg = "white"` to the Figure 2 save so the background is not transparent

## (d) Results and Figures
**Gene:** TPM2, "tropomyosin 2", Gene ID 7169, chromosome 9.

**Protein:** NP_003280.2, 284 amino acids. Residue 133 is Arg, consistent with p.Arg133Trp.

**Composition (top three):** Glu (60), Lys (39), Ala (34). Leu (29) is next. Pro and Trp are both absent (0).

| Figure | File |
|---|---|
| Figure 1 | `figures/01_amino_acid_composition.png` |
| Figure 2 | `figures/02_variant_position.png` |
| Figure 3 | `figures/protein_structure_or_splice_diagram.png` |

**Figure 1.** Amino-acid composition of NP_003280.2 (284 aa). Glu, Lys and Ala are most abundant.

**Figure 2.** Location of p.Arg133Trp (residue 133) on the 284-aa TPM2 reference protein (NP_003280.2). The marker shows sequence position only, not structural or functional effect.

**Figure 3.** AlphaFold DB predicted structure of TPM2 (UniProt P07951, AF-P07951-F1), accessed 2026-10-28, colored by pLDDT confidence. Residue 133 (Arg, site of p.Arg133Trp) is outlined in yellow. This is a single-chain computational prediction. It does not show the coiled-coil dimer, binding partners, or the effect of the Trp substitution. Source: AlphaFold DB (Jumper et al. 2021; Varadi et al. 2022).

## Step 7. Explain the molecular mechanism
1. Describe what your normal protein does in skeletal, cardiac, or smooth muscle (as appropriate). Identify whether it drives the actin–myosin cycle, regulates calcium, or maintains musclestructure.
- TPM2 encodes beta-tropomyosin, an actin-binding protein of skeletal muscle. Two chains form a coiled-coil dimer that lies along the actin thin filament. In relaxed muscle, tropomyosin blocks the myosin-binding sites on actin. When Ca2+ binds troponin C, the troponin complex shifts tropomyosin, which exposes the sites so myosin cross-bridges can form and the muscle contracts. TPM2 therefore regulates the actin-myosin cycle (thin-filament regulation).

2. Explain what is known about the assigned variant, separating observed experimental evidencefrom computational prediction and clinical association.
- See Discussion.

3. Write a short causal chain: gene → protein → molecular function →assigned variant →possibleeffect on contraction, relaxation, or muscle integrity.
- TPM2 gene -> beta-tropomyosin (coiled-coil dimer on actin) -> regulates myosin access to actin through the troponin-Ca2+ system -> c.397C>T replaces Arg133 with Trp -> altered tropomyosin behavior on the thin filament (supported by experiments, details uncertain) -> abnormal regulation of contraction and relaxation -> muscle weakness and joint contractures reported in distal arthrogryposis.

4. Propose one experiment that could test your interpretation, such as ATPase measurement, calcium imaging, binding assay, microscopy, or force measurement.
- See section (f).

### Discussion (Experimental vs. Clinical vs. Predicted)
TPM2 encodes beta-tropomyosin, which forms a coiled-coil dimer along the actin thin filament. The assigned variant replaces Arg133 with Trp. My composition count shows Trp is absent from this protein, so this change adds a residue type that is normally missing.

**Experimentally observed.** Ochala et al. (2007) tested single chemically skinned muscle fibres from two patients carrying R133W and from healthy controls. In type I fibres, the patient fibres produced less force for their size (lower specific force), recovered force more slowly after a release (lower ktr), and shortened faster when unloaded (higher V0). The force-pCa relationship was the same as in controls, so Ca2+ activation of the contractile proteins was not changed. The authors concluded that R133W slows cross-bridge attachment and speeds up detachment, so fewer myosin heads are in the strong actin-bound state, which fits muscle weakness without muscle wasting. This is evidence from patient fibres, not from a purified-protein assay. A ClinVar submitter comment (GeneDx) also cites a separate functional study (Ochala et al., 2010) reporting that the variant partly inhibits Ca2+- and myosin-induced tropomyosin movement over the thin filament; I did not read that paper. Marston et al. (2013) describe gain-of-function effects for TPM2 mutations in myopathy patients, so the effect may not be a simple loss of function.

**Clinical association.** ClinVar classifies the variant as Pathogenic. It was reported in families with distal arthrogryposis, including a mother and daughter (Tajsharghi et al. 2007) and a Korean family (Ko et al. 2013), and submitter comments describe de novo occurrences. The variant is absent from gnomAD according to submitters, which shows it is rare but does not prove causation. A ClinVar classification is a clinical interpretation based on evidence, not an experiment itself.

**Predicted.** Submitters report REVEL score 0.79 and 3Cnet score 0.93, and a modeling analysis predicting disrupted protein function. These are computational predictions. The AlphaFold model shows residue 133 inside a very-high-confidence helix (pLDDT above 90), but it is a single-chain model — so it cannot show the dimer interface, actin or troponin binding contacts, or how the Trp substitution alters those interactions.

**Uncertain.** The exact molecular mechanism is not settled. The data support altered thin-filament regulation, but whether Trp changes tropomyosin dimer stiffness, actin-binding affinity, or the magnitude and timing of the Ca²⁺ response remains unproven by direct measurement.

## (f) Limitations and Validation Experiment

**Limitations**
- Amino-acid composition cannot show structure, function, or variant effect.
- The AlphaFold model is a prediction of one chain; it does not show the dimer or binding partners, and it does not model the Trp substitution.
- ClinVar classifications and evidence can change; this record reflects status as of 2026-10-28.
- I did not independently verify the cited papers beyond the ClinVar citation list and submitter comments.

**Proposed experiment**
Measure force against Ca²⁺ concentration (pCa–force curves) in skinned muscle fibres reconstituted with recombinant wild-type or R133W beta-tropomyosin, comparing Ca²⁺ sensitivity (pCa₅₀), maximal force, and cooperativity. As a complementary test, an F-actin co-sedimentation assay or in vitro motility assay would measure actin binding and filament sliding.

**Questions to Answer**

| Question | Answer |
|---|---|
| What does your assigned gene encode, and how does this protein contribute to muscle function? | It encodes beta-tropomyosin, an actin-binding coiled-coil protein. In relaxed muscle it blocks myosin-binding sites on actin; when Ca²⁺ binds troponin C, tropomyosin shifts, exposing the sites so cross-bridges can form and the muscle can contract. |
| What was the exact protein accession and length retrieved with R? | NP_003280.2, 284 amino acids (NM_003289.4, Gene ID 7169). |
| What did the amino-acid composition chart show? What can it NOT tell you? | Glu (60), Lys (39) and Ala (34) were most abundant, with many charged residues — consistent with a helical coiled-coil. It cannot show structure, function, or the effect of the variant. |
| What is your fixed ClinVar variant, and what is its current classification and evidence strength? | c.397C>T (p.Arg133Trp), ClinVar 12463, **Pathogenic**, 9 of 12 submissions contributing, two-star review status (accessed 2026-10-28). |
| What is the scientifically supported mechanism linking this variant to muscle function, and what remains uncertain? | Experimental work supports altered regulation of contraction by R133W tropomyosin; the molecular details (loss vs. gain of function, effect on actin binding or Ca²⁺ response) remain uncertain. |
| What additional experiment would help test the proposed mechanism? | pCa–force measurements in skinned fibres with recombinant wild-type vs. R133W beta-tropomyosin, plus co-sedimentation or motility assays. |

**How to rerun** 
1. Install R and RStudio, then `install.packages(c("rentrez", "ggplot2"))`.
2. Open RStudio in the project folder and run `analysis.R` from top to bottom.
3. If an NCBI request fails, rerun only the failed step. Do not change the gene or variant.

## (g) Sources and Access Dates
All databases accessed: **2026-10-28**

| Source | Link |
|---|---|
| NCBI Gene — TPM2 (ID 7169) | https://www.ncbi.nlm.nih.gov/gene/7169 |
| NCBI Nucleotide — NM_003289.4 | https://www.ncbi.nlm.nih.gov/nuccore/NM_003289.4 |
| NCBI Protein — NP_003280.2 | https://www.ncbi.nlm.nih.gov/protein/NP_003280.2 |
| ClinVar — Variation ID 12463 | https://www.ncbi.nlm.nih.gov/clinvar/variation/12463 |
| UniProt — P07951 | https://www.uniprot.org/uniprotkb/P07951 |
| AlphaFold DB — AF-P07951-F1 | https://alphafold.ebi.ac.uk/search/text/P07951 |
| rentrez Documentation | https://docs.ropensci.org/rentrez/ |

### References
- Tajsharghi H, et al. Distal arthrogryposis and muscle weakness associated with a beta-tropomyosin mutation. *Neurology*. 2007. PMID 17339586.
- Ochala J, et al. Effects of a R133W beta-tropomyosin mutation on regulation of muscle contraction in single human muscle fibres. *J Physiol*. 2007. PMID 17430991.
- Ko JM, et al. First Korean family with a mutation in TPM2 associated with Sheldon-Hall syndrome. *J Korean Med Sci*. 2013. PMID 23678273.
- Marston S, et al. Mutations in repeating structural motifs of tropomyosin cause gain of function in skeletal muscle myopathy patients. *Hum Mol Genet*. 2013. PMID 23886664.
- Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021.
- Varadi M, et al. AlphaFold Protein Structure Database. *Nucleic Acids Res*. 2022.
