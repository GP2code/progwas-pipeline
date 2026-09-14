# GitHub Actions for Automated Docker Builds

This repository includes GitHub Actions workflows that automatically build and push Docker images to GitHub Container Registry (ghcr.io) when you create version tags.

## Setup Instructions

### 1. Enable GitHub Container Registry

No additional setup required! The workflow uses `GITHUB_TOKEN` which is automatically provided by GitHub Actions.

### 2. Set Package Visibility (Optional)

After the first build, you may want to make the image public:

1. Go to your GitHub profile → **Packages**
2. Find the `progwas-pipeline` package
3. Click **Package settings**
4. Scroll to **Danger Zone**
5. Click **Change visibility** → **Public**

### 3. Trigger a Build

The workflow automatically triggers when you push a version tag:

```bash
# Create and push a version tag
git tag -a v2.0.0 -m "Release version 2.0.0"
git push origin v2.0.0
```

This will:
- Build Docker image from `Dockerfile.ubuntu22`
- Tag with:
  - The version number (e.g., `2.0.0`)
  - Major.minor version (e.g., `2.0`)
  - Major version (e.g., `2`)
  - `latest` (if on main branch)
  - `slim`
- Push to GitHub Container Registry (ghcr.io)
- Build for multiple architectures (amd64, arm64)

### 4. Manual Trigger (Optional)

You can also manually trigger the workflow:

1. Go to **Actions** tab in GitHub
2. Select **Build and Push Docker Image**
3. Click **Run workflow**
4. Select branch and click **Run workflow**

## Image Tags

After pushing `v1.0.0`, the following tags will be available:

- `ghcr.io/gp2code/progwas-pipeline:1.0.0`
- `ghcr.io/gp2code/progwas-pipeline:1.0`
- `ghcr.io/gp2code/progwas-pipeline:1`
- `ghcr.io/gp2code/progwas-pipeline:latest`

## Using Versioned Images

In your profile configs, you can now reference specific versions:

```groovy
process {
  container = 'ghcr.io/gp2code/progwas-pipeline:1.0.0'  // Specific version
    // or
  container = 'ghcr.io/gp2code/progwas-pipeline:1'      // Major version
    // or
  container = 'ghcr.io/gp2code/progwas-pipeline:latest' // Latest release
}
```

## Pulling Images

Images are public (after you set visibility) and can be pulled without authentication:

```bash
docker pull ghcr.io/gp2code/progwas-pipeline:1.0.0
```

For Singularity (on HPC):

```bash
singularity pull docker://ghcr.io/gp2code/progwas-pipeline:1.0.0
```

## Workflow Features

- ✅ Multi-architecture builds (amd64, arm64)
- ✅ Layer caching for faster builds
- ✅ Automatic semantic versioning
- ✅ Integrated with GitHub (no separate service)
- ✅ Manual trigger option
- ✅ Free unlimited storage for public images
- ✅ Automatic linking to repository

## Troubleshooting

### Build Fails

Check the **Actions** tab in GitHub to see the build logs. Common issues:

- **Permission denied**: Workflow needs `packages: write` permission (already configured)
- **Authentication failed**: GitHub token should work automatically
- **Build errors**: Check Dockerfile.ubuntu22 syntax

### Tag Not Triggering

The workflow only triggers on tags matching `v*.*.*` pattern:
- ✅ `v1.0.0`, `v2.1.3`
- ❌ `v1`, `1.0.0`, `release-1.0`

Make sure to use the correct format:
```bash
git tag -a v2.0.0 -m "Description"
git push origin v2.0.0
```
