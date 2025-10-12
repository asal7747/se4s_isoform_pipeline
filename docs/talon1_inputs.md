## Inputs used (exact file list)

# External inputs (source-provenanced)

ENCODE long-read datasets (PacBio Iso-Seq; mouse)
- ENCFF003OWX  (C2C12/myotube; used to produce work/longread/bulk/ENCFF003OWX/labeled/clean_labeled.sam)  [ENCODE]
- ENCFF019HRC  (C2C12; used to produce work/longread/bulk/ENCFF019HRC/labeled/clean_labeled.sam)  [ENCODE]
- ENCFF669LWV  (C2C12; used to produce work/longread/bulk/ENCFF669LWV/labeled/clean_labeled.sam)  [ENCODE]
- ENCFF676BYQ  (C2C12; used to produce work/longread/bulk/ENCFF676BYQ/labeled/clean_labeled.sam)  [ENCODE]
Notes: These ENCODE accessions are PacBio Iso-Seq mouse data; the run used the four accessions above as the starting reads for labeling and subsequent TALON annotation.  [Standards: ENCODE]

Reference genome (mm10/GRCm38)
- mm10.fa.gz  (genome FASTA; UCSC bigZips) → refs/mm10.fa(.gz)
- mm10.fa.fai (samtools index generated locally)
Notes: mm10/GRCm38 was used as the labeling reference and read-wise coordinate frame.  [UCSC bigZips]

GENCODE annotation (vM21; mm10)
- gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz  → used to build and interpret mm10_vM21.db
Notes: vM21 is the GENCODE M21 release for mouse (mm10/GRCm38), matching the TALON database assembly and annotation version.  [GENCODE vM21]

TALON database (prebuilt)
- mm10_vM21.db  (TALON SQLite database used for annotate and exporter steps)
Notes: The annotation run updates this database and then writes the unified read-wise TSV.  [TALON]

Software/runtime
- TALON 5.0 (Bioconda) in an x86_64 (Rosetta) conda env with Python 3.7
- Docker Desktop (macOS) and image “kwellswrasman/talon” (optional, used for talon_label_reads)
Notes: Labeling was executed in the container; annotation/export used the local x86_64 conda environment.  [Bioconda TALON][Docker docs]

Reference and database
- Reference FASTA: /Users/asheralbrecht/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/refs/mm10.fa  [required for talon_label_reads]  [web:5]
- TALON database (mm10 vM21): /Users/asheralbrecht/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/work/longread/talon/mm10_vM21.db  [web:5]

Labeled SAM inputs (per dataset)
- ENCFF003OWX: /Users/asheralbrecht/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/work/longread/bulk/ENCFF003OWX/labeled/clean_labeled.sam  [web:5]
- ENCFF019HRC: /Users/asheralbrecht/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/work/longread/bulk/ENCFF019HRC/labeled/clean_labeled.sam  [web:5]
- ENCFF669LWV: /Users/asheralbrecht/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/work/longread/bulk/ENCFF669LWV/labeled/clean_labeled.sam  [web:5]
- ENCFF676BYQ: /Users/asheralbrecht/Library/CloudStorage/GoogleDrive-asal7747@colorado.edu/My Drive/Fall Semester 2025/SE4S/work/longread/bulk/ENCFF676BYQ/labeled/clean_labeled.sam  [web:5]

Optional base alignment inputs (only if re-labeling was needed)
- Per dataset BAM (one of): clean.bam or aligned.bam under /Users/.../work/longread/bulk/<DATASET>/  [used to derive clean.sam if missing]  [web:5]

Generated during run (for provenance)
- TALON config used for annotation: $HOME/Desktop/talon_scratch/config_local.csv  [four columns: dataset, description, platform, SAM path]  [web:5]
- TALON outprefix and log: $HOME/Desktop/talon_scratch/bulk_run_local_*  [read_annot TSV and QC log produced here]  [web:5]
