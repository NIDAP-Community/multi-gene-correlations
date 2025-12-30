# Quick Start Guide

## Push to GitHub Container Registry

### 1. Login to GHCR

```bash
# Set your GitHub username and token
export GITHUB_USER="your-github-username"
export GITHUB_TOKEN="your-personal-access-token"

# Login
echo $GITHUB_TOKEN | docker login ghcr.io -u $GITHUB_USER --password-stdin
```

### 2. The image is already built and tagged:

```bash
ghcr.io/nidap-community/multi-gene-correlations:latest
```

### 3. Push to Registry

```bash
docker push ghcr.io/nidap-community/multi-gene-correlations:latest
```

### 4. Optional: Add version tag

```bash
docker tag ghcr.io/nidap-community/multi-gene-correlations:latest \
  ghcr.io/nidap-community/multi-gene-correlations:1.0.0

docker push ghcr.io/nidap-community/multi-gene-correlations:1.0.0
```

### 5. Make Package Public

1. Go to https://github.com/orgs/NIDAP-Community/packages
2. Find `multi-gene-correlations`
3. Package settings → Change visibility → Public

## Test the Published Image

```bash
# Logout to test public access
docker logout ghcr.io

# Pull the image
docker pull ghcr.io/nidap-community/multi-gene-correlations:latest

# Run a test
docker run --rm ghcr.io/nidap-community/multi-gene-correlations:latest
```

## Create Singularity Image

```bash
singularity pull multi-gene-correlations.sif \
  docker://ghcr.io/nidap-community/multi-gene-correlations:latest
```

## Local Testing Commands

```bash
# Test with your data
docker run --rm \
  -v "/Users/maggiec/GitHub/Maggie/Harris/Data:/data" \
  -v "$PWD/output:/output" \
  ghcr.io/nidap-community/multi-gene-correlations:latest \
  --counts /data/IODC-0001-TCGA-PRAD.normalized_counts.tsv \
  --metadata /data/IODC-0001-TCGA-PRAD.meta.sample.csv \
  --gene_column gene \
  --sample_column sample2 \
  --genes "GATM,GAMT,CKM,SLC6A8,AR,PTEN,SPRY2" \
  --category_column samptype \
  --categories "TP,NT" \
  --signature_name "Creatine" \
  --signature_genes "PTEN,SPRY2,CKM" \
  --output_dir /output

# Shell access for debugging
docker run --rm -it \
  -v "/Users/maggiec/GitHub/Maggie/Harris/Data:/data" \
  --entrypoint bash \
  ghcr.io/nidap-community/multi-gene-correlations:latest
```

## View Documentation in Container

```bash
docker run --rm --entrypoint cat \
  ghcr.io/nidap-community/multi-gene-correlations:latest \
  /app/README.md
```
