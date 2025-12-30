# Building and Publishing to GitHub Container Registry

This guide explains how to build and publish this container to the NIDAP-Community GitHub Container Registry.

## Prerequisites

1. **GitHub Account** with access to the NIDAP-Community organization
2. **Docker** installed locally
3. **GitHub Personal Access Token (PAT)** with `write:packages` and `read:packages` permissions

### Creating a GitHub PAT

1. Go to GitHub Settings → Developer settings → Personal access tokens → Tokens (classic)
2. Click "Generate new token (classic)"
3. Select scopes:
   - `write:packages` (to publish)
   - `read:packages` (to pull)
   - `delete:packages` (optional, for cleanup)
4. Save the token securely

## Step 1: Login to GitHub Container Registry

```bash
echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin
```

Replace `USERNAME` with your GitHub username and set `GITHUB_TOKEN` environment variable to your PAT.

## Step 2: Build the Image

```bash
cd /Users/maggiec/GitHub/Maggie/NIDAP/Templates/Correlation_Plot

docker build -t ghcr.io/nidap-community/multi-gene-correlations:latest .
```

Optional: Tag with version number:
```bash
docker tag ghcr.io/nidap-community/multi-gene-correlations:latest \
  ghcr.io/nidap-community/multi-gene-correlations:1.0.0
```

## Step 3: Test the Image Locally

```bash
docker run --rm ghcr.io/nidap-community/multi-gene-correlations:latest --help
```

Run a full test with your data:
```bash
docker run --rm \
  -v "/path/to/data:/data" \
  -v "$PWD/output:/output" \
  ghcr.io/nidap-community/multi-gene-correlations:latest \
  --counts /data/counts.tsv \
  --metadata /data/metadata.csv \
  --gene_column Gene \
  --sample_column SampleID \
  --genes "PTEN,CKM" \
  --category_column Category \
  --categories "A,B" \
  --signature_genes "PTEN,CKM" \
  --output_dir /output
```

## Step 4: Push to GitHub Container Registry

Push the latest tag:
```bash
docker push ghcr.io/nidap-community/multi-gene-correlations:latest
```

Push versioned tag (if created):
```bash
docker push ghcr.io/nidap-community/multi-gene-correlations:1.0.0
```

## Step 5: Make Package Public

1. Go to https://github.com/orgs/NIDAP-Community/packages
2. Find `multi-gene-correlations`
3. Click on the package
4. Go to "Package settings"
5. Scroll to "Danger Zone" → "Change visibility"
6. Select "Public"

## Step 6: Update Repository Settings

In the GitHub repository:

1. Go to the repository on GitHub
2. Click "Settings" → "Actions" → "General"
3. Under "Workflow permissions", ensure "Read and write permissions" is selected

## Step 7: Link Package to Repository

1. Go to the package page: https://github.com/orgs/NIDAP-Community/packages/container/multi-gene-correlations
2. Click "Connect repository"
3. Select the `multi-gene-correlations` repository

## Verification

Test that others can pull the image:

```bash
docker logout ghcr.io
docker pull ghcr.io/nidap-community/multi-gene-correlations:latest
```

## Building Singularity Image

After publishing to ghcr.io, users can build Singularity images:

```bash
singularity pull multi-gene-correlations.sif \
  docker://ghcr.io/nidap-community/multi-gene-correlations:latest
```

Or build from local Docker:
```bash
singularity build multi-gene-correlations.sif \
  docker-daemon://ghcr.io/nidap-community/multi-gene-correlations:latest
```

## Automated Publishing with GitHub Actions

Create `.github/workflows/docker-publish.yml` in your repository:

```yaml
name: Publish Docker Image

on:
  push:
    branches: [ main ]
    tags: [ 'v*.*.*' ]
  pull_request:
    branches: [ main ]

env:
  REGISTRY: ghcr.io
  IMAGE_NAME: ${{ github.repository }}

jobs:
  build-and-push:
    runs-on: ubuntu-latest
    permissions:
      contents: read
      packages: write

    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Log in to Container Registry
        uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ${{ env.REGISTRY }}/${{ env.IMAGE_NAME }}

      - name: Build and push Docker image
        uses: docker/build-push-action@v5
        with:
          context: .
          push: ${{ github.event_name != 'pull_request' }}
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
```

This will automatically build and push images on every commit to main and on version tags.

## Troubleshooting

### Permission Denied

Ensure your PAT has the correct scopes and you're logged in:
```bash
docker login ghcr.io
```

### Package Already Exists

If the package exists but you can't push:
1. Check organization membership
2. Verify package permissions
3. Contact organization admin

### Image Too Large

Check image size:
```bash
docker images ghcr.io/nidap-community/multi-gene-correlations
```

Optimize if needed by removing unnecessary files or using multi-stage builds.

## Support

For issues with publishing, contact the NIDAP-Community administrators or open an issue in the repository.
