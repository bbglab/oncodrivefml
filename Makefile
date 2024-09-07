ROOT_DIR := $(shell echo $(dir $(lastword $(MAKEFILE_LIST))) | sed 's|/*$$||')

SHELL := /bin/bash

DOCS_DIR := $(ROOT_DIR)/docs
VENV_DIR := $(ROOT_DIR)/.venv

define version
$(shell uv run python -c "from oncodrivefml import __version__; print(__version__)")
endef

define git_tag_or_sha
$(shell git describe --tags --exact-match 2>/dev/null || git rev-parse --short HEAD)
endef

define image
bbglab/oncodrivefml:$(call version)
endef

IMAGE_FILE := oncodrivefml.tar

BOLDRED := $(shell tput bold && tput setaf 1)
BOLDGREEN := $(shell tput bold && tput setaf 2)
BOLDYELLOW := $(shell tput bold && tput setaf 3)
BOLDBLUE := $(shell tput bold && tput setaf 4)
LIGHTBLUE := $(shell tput setaf 6)
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
	@echo "$(BOLDGREEN)  build-image  $(WHITE)-> Build the Docker image"
	@echo "$(BOLDGREEN)  docker-login $(WHITE)-> Log in to DockerHub"
	@echo "$(BOLDGREEN)  push-image   $(WHITE)-> Push the Docker image into DockerHub"
	@echo "$(BOLDGREEN)  docs         $(WHITE)-> Generate the documentation"
	@echo "$(BOLDGREEN)  run-example  $(WHITE)-> Run the included example using the Docker image"
	@echo "$(BOLDGREEN)  clean        $(WHITE)-> Clean the working directory (build files, virtual environments, caches)"
	@echo "$(RESET)"

.PHONY: uv-installed
uv-installed:
	@if ! which uv > /dev/null; then \
		echo "$(BOLDRED)This project build is managed by $(BOLDYELLOW)uv$(BOLDRED), which is not installed.$(RESET)"; \
		echo "$(LIGHTBLUE)Please follow these instructions to install it:$(RESET)"; \
		echo "$(LIGHTBLUE)--> $(BOLDBLUE)https://docs.astral.sh/uv/#getting-started$(RESET)"; \
		exit 1; \
	fi

.PHONY: ruff-installed
ruff-installed: uv-installed
	@if ! which ruff > /dev/null; then \
		echo "$(BOLDRED)This project requires $(BOLDYELLOW)ruff$(BOLDRED), which is not installed.$(RESET)"; \
		echo "$(LIGHTBLUE)Installing it with $(BOLDYELLOW)uv tool install ruff$(RESET)"; \
		uv tool install ruff; \
		ruff --version; \
	fi

.PHONY: sphinx-installed
sphinx-installed: uv-installed
	@if ! uv pip show sphinx > /dev/null; then \
		echo "$(BOLDRED)This project requires $(BOLDYELLOW)sphinx$(BOLDRED), which is not installed.$(RESET)"; \
		echo "$(LIGHTBLUE)Installing it with $(BOLDYELLOW)uv pip install --requirements optional-requirements.txt$(RESET)"; \
		uv venv; \
		uv pip install --requirements optional-requirements.txt; \
		uv run sphinx-build --version; \
	fi

.PHONY: checks
checks: check-format check-lint check-docker

.PHONY: check-format
check-format: ruff-installed
	@echo "$(BOLDGREEN)Checking code format ...$(RESET)"
	ruff format --check
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: check-lint
check-lint: ruff-installed
	@echo "$(BOLDGREEN)Checking lint ...$(RESET)"
	ruff check
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
check-version: uv-installed
	@echo "$(BOLDGREEN)Checking that the version matches the tag ...$(RESET)"
	@if [ "$(call version)" != "$(call git_tag_or_sha)" ]; then \
	    echo "$(BOLDRED)==> Version $(BOLDYELLOW)$(call version)$(BOLDRED) doesn't match the git tag $(BOLDYELLOW)$(call git_tag_or_sha)$(BOLDRED) !!!$(RESET)"; \
		echo "$(BOLDRED)==> Please update the $(BOLDYELLOW)__version__$(BOLDRED) in $(BOLDYELLOW)oncodrivefml/__init__.py$(BOLDRED) and re-create the tag.$(RESET)"; \
	    exit 1; \
	fi
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: format
format: ruff-installed
	@echo "$(BOLDGREEN)Formatting code ...$(RESET)"
	ruff format

.PHONY: build-dist
build-dist: uv-installed
	@echo "$(BOLDGREEN)Building packages ...$(RESET)"
	uv build

.PHONY: build-sdist
build-sdist:
	@echo "$(BOLDGREEN)Building sdist package ...$(RESET)"
	uv build --sdist

.PHONY: build-wheels
build-wheels: uv-installed
	@echo "$(BOLDGREEN)Building wheels ...$(RESET)"
	uv venv
	uv pip install cibuildwheel
	uv run cibuildwheel --output-dir dist

.PHONY: publish-dist
publish-dist: uv-installed
	@echo "$(BOLDGREEN)Publishing OncodriveCLUSTL $(BOLDYELLOW)$(call version)$(BOLDGREEN) to PyPI ...$(RESET)"
	@if [ -z "$(PYPI_TOKEN)" ]; then \
		echo "$(BOLDRED)==> Missing PyPI token !!!$(RESET)"; \
		exit 1; \
	fi
	uv publish --token $(PYPI_TOKEN) dist/**/*

.PHONY: docker-login
docker-login:
	@echo "$(BOLDGREEN)Log in to DockerHub ...$(RESET)"
	@if [[ -z "$(DOCKER_USERNAME)" || -z "$(DOCKER_PASSWORD)" ]]; then \
		echo "$(BOLDRED)==> Missing DockerHub credentials !!!$(RESET)"; \
		exit 1; \
	fi
	@(echo "$(DOCKER_PASSWORD)" | docker login -u $(DOCKER_USERNAME) --password-stdin) || (echo "$(BOLDRED)==> Failed to log in !!!$(RESET)"; exit 1)

.PHONY: build-image
build-image: uv-installed
	@echo "$(BOLDGREEN)Building Docker image $(BOLDYELLOW)$(call image)$(BOLDGREEN) ...$(RESET)"
	docker build --progress=plain -t $(call image) .
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: save-image
save-image: build-image
	@echo "$(BOLDGREEN)Saving Docker image $(BOLDYELLOW)$(call image)$(BOLDGREEN) ...$(RESET)"
	docker save -o $(IMAGE_FILE) $(call image)
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: load-image
load-image: uv-installed
	@echo "$(BOLDGREEN)Loading Docker image $(BOLDYELLOW)$(call image)$(BOLDGREEN) ...$(RESET)"
	docker load -i $(IMAGE_FILE)
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: build-image
push-image: uv-installed
	@echo "$(BOLDGREEN)Pushing the Docker image into the DockerHub ...$(RESET)"
	docker push $(call image)
	@echo "$(BOLDGREEN)==> Success!$(RESET)"

.PHONY: docs
docs: sphinx-installed
	(source .venv/bin/activate; make -C $(DOCS_DIR) html)

.PHONY: run-example
run-example: build-image
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
	rm -rf ./oncodrivefml.egg-info ./dist ./.ruff_cache ./.eggs $(DOCS_DIR)/build ./.venv
	find oncodrivefml \( -name '*.c' -o -name '*.so' \) -type f -exec rm {} +
	find . -name "__pycache__" -type d -exec rm -r {} +
