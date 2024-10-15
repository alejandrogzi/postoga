#!/bin/bash

# author = "Alejandro Gonzales-Irribarren"
# email = "alejandrxgzi@gmail.com"
# github = "https://github.com/alejandrogzi"
# version: 0.9.3-devel


TEST_PY="postoga_test.py"
TEST_PROJECT="test_project"
ENV_YML="env.yml"
PISCO_ENV="postoga"
CONDA_FLAG=false

echo ""
echo ">> configuring environment to run postoga"

for arg in "$@"
do
    if [ "$arg" == "--conda" ]; then
        CONDA_FLAG=true
        break
    fi
done


if $CONDA_FLAG; then
  echo ">> using conda to create environment..."

  if ! command -v conda &> /dev/null
  then
      echo "!ERROR: conda could not be found, please install Miniconda or Anaconda."
      exit
  fi

  if [ -f $ENV_YML ]; then
      echo ">> creating conda environment..." && echo ""
      conda env create -f $ENV_YML

      if [[ $? -eq 0 ]]; then
          sleep 3
          eval "$(conda shell.bash hook)"
          conda activate $PISCO_ENV
      else
          echo ">> conda env create failed, activating existing environment..."
          sleep 3
          eval "$(conda shell.bash hook)"
          conda activate $PISCO_ENV
      fi
  else
      echo "!ERROR: $ENV_YML not found. Please clone the repository again!"
      exit
  fi
else
  echo ">> using hatch to create environment..."

  hatch env create

  if [[ $? -ne 0 ]]; then
      echo "!ERROR: hatch env create failed. Please check the error message above."
      exit 1
  fi

  VENV_PATH=$(hatch env find)
  if [[ $? -ne 0 ]]; then
      echo "!ERROR: Failed to find the hatch environment."
      exit 1
  fi

  source "$VENV_PATH/bin/activate"
  if [[ $? -ne 0 ]]; then
      echo "!ERROR: Failed to activate the hatch environment."
      exit 1
  fi

  echo ">> environment created successfully!"
fi

echo ">> building rust extensions and binaries..."
cd rustools && maturin develop --release
cd ..
cargo install bed2gtf --force
cargo install bed2gff --force
cargo install gxf2bed --force
echo ">> rust extensions and binaries built successfully!"

echo "" && echo ">> do you want to run the test? [y/n]"
read -r run_test

if [ "$run_test" == "y" ]; then
    echo ">> running test..."
    make test

    if [[ $? -ne 0 ]]; then
        echo "!ERROR: test failed. Please check the error message above."
        exit
    fi

    # clean up
    rm -f $TEST_PY && rm -rf $TEST_PROJECT
fi

echo "" && echo ">> SUCCESS: environment has been configured successfully. You are now ready to use postoga!"
