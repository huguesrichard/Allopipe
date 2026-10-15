# Running the tests

Run the commands below from the repository root. The suite tests the Python sources
and Nextflow workflow in the current checkout.

## Environment setup

Keep `allopipe` for application runs and `allopipe-test` for tests. The test
environment installs the runtime dependencies plus the test tools separately,
using two dependency files:

- [`allopipe.yml`](../recipes/conda/allopipe.yml) defines the runtime dependencies,
  shared with the Nextflow pipeline.
- [`dev-packages.txt`](../recipes/conda/dev-packages.txt) adds pytest and pytest-cov
  using Conda package specifications.

For a fresh setup, create and activate the environment, then add the test tools:

```sh
conda env create -n allopipe-test -f recipes/conda/allopipe.yml
conda activate allopipe-test
conda install --override-channels -c conda-forge \
  --freeze-installed --file recipes/conda/dev-packages.txt
```

If `allopipe-test` already exists, activate it with `conda activate allopipe-test`.
Install the test tools only if they have not been added yet; otherwise, go directly
to running the suite.

`--override-channels -c conda-forge` restricts this installation to conda-forge,
matching the intent of `nodefaults` in the runtime YAML. `--freeze-installed` keeps
installed dependencies unchanged while allowing additional dependencies required
by the test tools. Add `--dry-run` to preview the installation.

A previous Nextflow run alone does not create either named environment. Nextflow
manages its own cached Conda environments from `allopipe.yml` and does not read
`dev-packages.txt`. Create `allopipe-test` explicitly for testing, even if the
pipeline has already run. The dependency files constrain package versions but do
not provide a complete environment lock.

## Run the suite

The pipeline prerequisites must already be installed (see the
[installation instructions](../README.md#requirements)). The Nextflow tests run
offline using a copy of the runtime cache. If Nextflow has never been launched,
run `nextflow -version` once to initialize that cache.

With the test environment activated:

```sh
python -m pytest -ra --durations=10
```

`-ra` lists non-passing outcomes, including skipped tests; `--durations=10` reports
the ten slowest tests. To run a single test file, append its path:

```sh
python -m pytest tests/test_parsing_functions.py -ra
```

## Coverage

The default test command does not collect coverage or update `htmlcov`. To generate
both a terminal summary and an HTML report:

```sh
python -m pytest -ra --cov=src --cov-report=term-missing --cov-report=html
```

Open `htmlcov/index.html` in a browser. This report measures Python code executed
under coverage; it does not measure Nextflow execution or automatically include
Python subprocesses launched by the workflow.
