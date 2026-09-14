# Verily Workbench Setup

[Official Instructions](https://support.workbench.verily.com/docs/guides/workflows/nextflow/)

**Note for GP2 Workbench**:  
Although you can install Workbench CLI on your local machine, you can only access data through a running VM (app) in the same workspace. All commands should be run on the VM for GP2 analyses.

---

## Preparation steps on VWB

- Create a workspace on [Verily Workbench](https://workbench.verily.com/)
- Add resources to the workspace:
  - `+ New resource` → GCS Bucket for Nextflow run logs and output
    - Resource type: **GCS Bucket**
    - Resource ID: e.g., `nf_files`
    - Description: "Bucket for Nextflow run logs and output"
  - `+ Data from catalog` → GP2 tier 2 data (if needed)
- Create App instance to run Nextflow:
  - `+ Add repository` → progwas-pipeline repository (`https://github.com/GP2code/progwas-pipeline.git`)
  - `+ New app instance` → Jupyter Lab (default settings)

## Running on VWB

- Open the created Jupyter Lab app instance
- Open a terminal in Jupyter Lab
- The VM terminal is already authenticated with your Google account and has access to the workspace data
- VM mounts the Resources to but nextflow workflow require the full path not the mounted path. For example, if you created a GCS bucket resource with ID `nf_files`, the mounted path is `/mnt/nf_files`, but the full path is `gs://<your-bucket-name>`. Use the full path in your Nextflow parameters. 

**Check pre-set environment variables:**
```bash
echo $GOOGLE_CLOUD_PROJECT
echo $GOOGLE_SERVICE_ACCOUNT_EMAIL
echo $WORKBENCH_USER_EMAIL
```

**Set environment variables for the pipeline:**
```bash
export NXF_SYNTAX_PARSER=v1
export STORE_ROOT='gs://<your-bucket-name>'  # Bucket you created above
export PROJECT_NAME='testrun'                # Any name for your project folder
export TOWER_ACCESS_TOKEN='<your-token>'    # Get from https://cloud.seqera.io/tokens
export REFERENCE_DIR='gs://<your-bucket-name>/References' # Your Reference Data Locations
```

**Pro tip**: Save these to `~/.vwb.env` for reuse:
```bash
cat > ~/.vwb.env << 'EOF'
export NXF_SYNTAX_PARSER=v1
export STORE_ROOT='gs://your-bucket-name'
export PROJECT_NAME='my_analysis'
export TOWER_ACCESS_TOKEN='your-token'
export REFERENCE_DIR='gs://your-bucket-name/References'
EOF

# Load before each run
source ~/.vwb.env
```

**Run the pipeline:**
```bash
cd ~/repos/progwas-pipeline
git pull origin main  # Update to latest code

wb nextflow -log z_$(date +%Y%m%d_%H%M%S).log run main.nf \
  -profile gcb_vwb \
  -params-file ./example/params_YML/test_cs_linear.yml \
  -with-tower \
  -resume
```

**Command breakdown:**
- `wb nextflow` - Workbench-wrapped Nextflow command
- `-log z_$(date +%Y%m%d_%H%M%S).log` - Timestamped log file (sorted last in directory)
- `run main.nf` - Main pipeline script
- `-profile gcb_vwb` - Recommended Google Cloud Batch profile for Verily Workbench
- `-params-file params.yml` - Parameter configuration file
- `-with-tower` - Enable Seqera Platform monitoring
- `-resume` - Resume from last successful checkpoint

**Note**: Using `-with-tower` enables web monitoring at https://cloud.seqera.io, which is helpful for detailed logs and debugging when Cloud Logging permissions are limited.


**Monitor batch jobs:**
```bash
wb gcloud batch jobs list
```

---

## CLI-based setup (alternative)

Install the CLI: [Instructions](https://support.workbench.verily.com/docs/guides/cli/cli_install_and_run/)

### Authenticate and create workspace
```bash
wb auth login
wb workspace set --id=<workspace-id>
```

### Create storage for Nextflow outputs
```bash
wb resource create gcs-bucket --id=nf_files \
  --description="Bucket for Nextflow run logs and output"
```

### Add the progwas-pipeline repository
```bash
wb resource add-ref git-repo \
  --id=progwas-pipeline \
  --repo-url=https://github.com/GP2code/progwas-pipeline.git

# Check resources
wb resource list
```

### Create JupyterLab app
```bash
wb app create gcp --app-config=jupyter-lab \
  --id=nextflow-jupyterlab \
  --description="JupyterLab notebook for running Nextflow"
```

### Run
Open JupyterLab via the web UI and follow the "Running on VWB" section above.

