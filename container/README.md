# Container Build and Runtime

This directory defines how the tool is containerized and executed.

This tool follows the **Container-as-a-Function (CAF)** model:
- non-interactive execution
- explicit inputs (CLI flags and mounted files)
- explicit outputs written to a user-specified directory
- no hidden state

The container is a runtime artifact, not a development shell.

---

## Canonical entrypoint for this repository

This repository’s canonical entrypoint is:

- `entrypoint.R` (repository root)

Core computation is implemented in:

- `Multi_Gene_Correlations_to_Signature.R` (repository root)

Note: This repository predates the newer `scripts/run_all.R` convention used in newer tool repos.
The CAF contract still applies exactly the same way.

---

## Canonical distribution model

This repository does NOT store container binaries in git.

### Canonical artifact
- A Docker/OCI image published to GitHub Container Registry (GHCR)

### Supported HPC runtime
- Singularity/Apptainer is supported by pulling and converting the Docker image:

singularity pull multi-gene-correlations.sif docker://ghcr.io/<org>/<tool>:<tag>


### Optional release artifact
- A prebuilt `.sif` may be attached to a GitHub Release for offline / restricted-network HPC use.
- `.sif` files are never committed to the repository.

---

## Platform assumptions

- Target runtime (Biowulf): **linux/amd64 (x86_64)**
- Canonical build host: **x86_64 VM with Docker**
- Apple Silicon (ARM64):
- supported for development
- container builds must explicitly target `linux/amd64`
- not the canonical build path

---

## Build options

### Option A — Docker-available environment (recommended)

docker build -t <org>/<tool>:dev -f container/Dockerfile .

During the build the Dockerfile preinstalls the tidyverse runtime stack
(`readr`, `forcats`, `vroom`, `hms`, `tzdb`, `bit`, `bit64`, `progress`,
`prettyunits`, `cpp11`, etc.), copies the synthetic data under `test_data/`,
and runs [`scripts/generate_clinical_correlation_plots.R`](../scripts/generate_clinical_correlation_plots.R)
as a smoke test to ensure all dependencies are satisfied before publishing.

Publish/tag according to project release conventions (e.g., GHCR).

---

### Option B — Apple Silicon (optional, non-canonical)

If building on Apple Silicon, explicitly target amd64:
docker buildx build
--platform linux/amd64
-t <org>/<tool>:dev
-f container/Dockerfile .


This may be slower due to emulation.

---

### Option C — HPC-restricted environment (Biowulf)

Biowulf does NOT support privileged container builds:
- no root
- no fakeroot
- limited unprivileged `proot` builds for simple definition files only

Standard approach:
- Build off-cluster (VM or local machine)
- Generate a `.sif` off-cluster
- Transfer `.sif` to Biowulf
- Run jobs using the `.sif`

Building containers directly on Biowulf is not the default workflow.

---

## Running the container

Run the tool using the container runtime you have available. The container must be invoked
non-interactively with explicit inputs and an explicit output directory.

### Docker (example)

docker run --rm
-v "<host_input_dir>:/data:ro"
-v "<host_output_dir>:/output"
ghcr.io/<org>/<tool>:<tag>
Rscript /work/entrypoint.R
--input_dir /data
--output_dir /output


### Singularity/Apptainer (example)
singularity pull <tool>.sif docker://ghcr.io/<org>/<tool>:<tag>

singularity run
--bind <host_input_dir>:/data
--bind <host_output_dir>:/output
<tool>.sif
--input_dir /data
--output_dir /output


(Exact flags depend on the entrypoint’s CLI; consult the top-level README for tool usage.)

### Helper wrapper for Singularity

[run_clinical_plots.sh](run_clinical_plots.sh) standardizes Singularity invocations by
binding the host data/output directories to `/data` and `/output`, setting
`R_LIBS_USER=/usr/local/lib/R/site-library`, and preserving CAF semantics. Example:

```
DATA_DIR=$PWD/test_data OUTPUT_DIR=$PWD/test_output \
container/run_clinical_plots.sh <tool>.sif --counts /data/counts.tsv ...
```

Customize the `DATA_DIR`, `OUTPUT_DIR`, and `EXTRA_BINDS` environment variables to
match the host filesystem layout.



