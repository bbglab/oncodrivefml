ROOT_DIR := $(shell echo $(dir $(lastword $(MAKEFILE_LIST))) | sed 's|/*$$||')

SHELL := /bin/bash

DOCS_DIR := $(ROOT_DIR)/docs
VENV_DIR := $(ROOT_DIR)/.venv

define version
$(shell $(VENV_DIR)/bin/python -c 'from oncodrivefml import __version__; print(__version__)')
endef

define image
bbglab/oncodrivefml:$(version)
endef

GIT_TAG_OR_SHA = $(shell git describe --tags --exact-match 2>/dev/null || git rev-parse --short HEAD)

BOLDRED := $(shell tput bold && tput setaf 1)
BOLDGREEN := $(shell tput bold && tput setaf 2)
BOLDYELLOW := $(shell tput bold && tput setaf 3)
WHITE := $(shell tput sgr0 && tput setaf 7)
RESET := $(shell tput sgr0)


.PHONY: help
help:
	@echo "$(BOLDYELLOW)Available targets:$(RESET)"
	@echo
	@echo "$(BOLDGREEN)  checks       $(WHITE)-> Run all the checks (format and lint)"
	@echo "$(BOLDGREEN)  check-format $(WHITE)-> Check for formatting errors"
	@echo "$(BOLDGREEN)  check-lint   $(WHITE)-> Check for lint errors"
	@echo "$(BOLDGREEN)  check-docker $(WHITE)-> Check the Dockerfile"
	@echo "$(BOLDGREEN)  format       $(WHITE)-> Format source code"
	@echo "$(BOLDGREEN)  build-dist   $(WHITE)-> Build source and wheel distribution files"
	@echo "$(BOLDGREEN)  install-dev  $(WHITE)-> Install the packages in editable mode"
	@echo "$(BOLDGREEN)  create-env   $(WHITE)-> Create a virtual environment"
	@echo "$(BOLDGREEN)  remove-env   $(WHITE)-> Remove the virtual environment"
	@echo "$(BOLDGREEN)  build-image  $(WHITE)-> Build the Docker image"
	@echo "$(BOLDGREEN)  docker-login $(WHITE)-> Log in to DockerHub"
	@echo "$(BOLDGREEN)  push-image   $(WHITE)-> Push the Docker image into DockerHub"
	@echo "$(BOLDGREEN)  docs         $(WHITE)-> Generate the documentation"
	@echo "$(BOLDGREEN)  run-example  $(WHITE)-> Run the included example using the Docker image"
	@echo "$(BOLDGREEN)  clean        $(WHITE)-> Clean the working directory (build files, virtual environments, caches)"
	@echo "$(RESET)"

$(VENV_DIR):
	@echo "$(BOLDYELLOW)Preparing virtual environment ...$(RESET)"
	python -m venv $(VENV_DIR)
	$(VENV_DIR)/bin/pip install -U pip ruff setuptools wheel build twine
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: checks
checks: check-format check-lint check-docker

.PHONY: check-format
check-format: $(VENV_DIR)
	@echo "$(BOLDGREEN)Checking code format ...$(RESET)"
	$(VENV_DIR)/bin/ruff format --check
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: check-lint
check-lint: $(VENV_DIR)
	@echo "$(BOLDGREEN)Checking lint ...$(RESET)"
	$(VENV_DIR)/bin/ruff check
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: check-docker
check-docker:
	@echo "$(BOLDGREEN)Checking Dockerfile ...$(RESET)"
	docker run --rm -i \
		-v $$(pwd):/project \
		hadolint/hadolint hadolint \
		--config /project/.hadolint.yaml \
		/project/Dockerfile
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: check-version
check-version: $(VENV_DIR)
	@echo "$(BOLDGREEN)Checking that the version matches the tag ...$(RESET)"
	@if [ "$(version)" != "$(GIT_TAG_OR_SHA)" ]; then \
	    echo "$(BOLDRED)==> Version $(BOLDYELLOW)$(version)$(BOLDRED) doesn't match the git tag $(BOLDYELLOW)$(GIT_TAG_OR_SHA)$(BOLDRED) !!!$(RESET)"; \
		echo "$(BOLDRED)==> Please update the $(BOLDYELLOW)__version__$(BOLDRED) in $(BOLDYELLOW)oncodrivefml/__init__.py$(BOLDRED) and re-create the tag.$(RESET)"; \
	    exit 1; \
	fi
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: format
format: $(VENV_DIR)
	@echo "$(BOLDGREEN)Formatting code ...$(RESET)"
	$(VENV_DIR)/bin/ruff format

.PHONY: build-dist
build-dist: $(VENV_DIR)
	@echo "$(BOLDGREEN)Building packages ...$(RESET)"
	$(VENV_DIR)/bin/python -m build

.PHONY: publish-dist
publish-dist: $(VENV_DIR)
	@echo "$(BOLDGREEN)Publishing OncodriveFML $(BOLDYELLOW)$(version)$(BOLDGREEN) to PyPI ...$(RESET)"
	@if [[ -z "$(PYPI_TOKEN)" ]]; then \
		echo "$(BOLDRED)==> Missing PYPI_TOKEN !!!$(RESET)"; \
		exit 1; \
	fi
	@$(VENV_DIR)/bin/twine upload --username __token__ --password $(PYPI_TOKEN) dist/*

.PHONY: install
install-dev: $(VENV_DIR)
	$(VENV_DIR)/bin/pip install -e .

.PHONY: create-env
create-env: $(VENV_DIR)

.PHONY: remove-env
remove-env:
	@echo "$(BOLDGREEN)Removing virtual environment ...$(RESET)"
	rm -rf $(VENV_DIR)

.PHONY: build-image
build-image: $(VENV_DIR)
	@echo "$(BOLDGREEN)Building Docker image $(BOLDYELLOW)$(image)$(BOLDGREEN) ...$(RESET)"
	docker build --progress=plain -t $(image) .
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: docker-login
docker-login:
	@echo "$(BOLDGREEN)Log in to DockerHub ...$(RESET)"
	@if [[ -z "$(DOCKER_USERNAME)" || -z "$(DOCKER_PASSWORD)" ]]; then \
		echo "$(BOLDRED)==> Missing DockerHub credentials !!!$(RESET)"; \
		exit 1; \
	fi
	@(echo "$(DOCKER_PASSWORD)" | docker login -u $(DOCKER_USERNAME) --password-stdin) || (echo "$(BOLDRED)==> Failed to log in !!!$(RESET)"; exit 1)

.PHONY: build-image
push-image: $(VENV_DIR)
	@echo "$(BOLDGREEN)Pushing the Docker image into the DockerHub ...$(RESET)"
	docker push $(image)
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: docs
docs: $(VENV_DIR)
	@if ! which $(VENV_DIR)/bin/sphinx-build > /dev/null; then \
		$(VENV_DIR)/bin/pip install -r optional-requirements.txt; \
	fi
	(source $(VENV_DIR)/bin/activate; make -C $(DOCS_DIR) html)

.PHONY: run-example
run-example:
	@echo "$(BOLDGREEN)Running example ...$(RESET)"
	docker run --rm -i \
		-v $${BGDATA_LOCAL:-$${HOME}/.bgdata}:/root/.bgdata \
		-v $$(pwd)/example:/data \
		--workdir /data \
		$(image) -i paad.txt.gz -e cds.tsv.gz --signature-correction wx --seed 123 --force
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: clean
clean:
	@echo "$(BOLDGREEN)Cleaning the repository ...$(RESET)"
	rm -rf ./oncodrivefml.egg-info ./dist $(VENV_DIR) ./.ruff_cache ./.eggs $(DOCS_DIR)/build
	find oncodrivefml \( -name '*.c' -o -name '*.so' \) -type f -exec rm {} +
	find . -name "__pycache__" -type d -exec rm -r {} +
