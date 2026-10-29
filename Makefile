.PHONY: help test fixture run run-demo validate clean

SERVICE_NAME ?= processor-seizure-ei
PYTHON ?= python3
IMAGE ?= processor-seizure-ei:dev
DEMO ?= ../workflow-epilepsy-ieeg/data/work/prep/ieeg
ONSETS ?=
TRUTH ?= $(DEMO)
OUT ?= data/output-demo
ANNOTATIONS ?=

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test      - pytest: synthetic BIDS-iEEG fixture through the processor (no Docker)"
	@echo "make fixture   - write the synthetic fixture to data/fixture/ (linked/ieeg = an ieeg-prep derivative)"
	@echo "make run       - build the image and run it in Docker over data/fixture/linked (ieeg/ subfolder), writing to data/output/"
	@echo "make run-demo  - run the image over the demo data (DEMO=$(DEMO), optional ONSETS=<seizure-dynasd output>) into OUT=$(OUT)"
	@echo "make validate  - score OUT against the expert labels (TRUTH=$(TRUTH))"
	@echo "make clean     - remove data/output*/ and data/fixture/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	$(PYTHON) -m tests.make_fixture data/fixture

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

run-demo:
	docker build --platform linux/amd64 -t $(IMAGE) .
	rm -rf "$(OUT)" && mkdir -p "$(OUT)"
	docker run --rm --platform linux/amd64 --cpus 1 --memory 2g --env-file dev.env \
		-v "$(abspath $(DEMO))":/data/input/ieeg:ro $(if $(ONSETS),-v "$(abspath $(ONSETS))":/data/input/onsets:ro) \
		-v "$(abspath $(OUT))":/data/output $(IMAGE)

validate:
	$(PYTHON) scripts/validate.py "$(OUT)" "$(abspath $(TRUTH))" $(ANNOTATIONS)

clean:
	rm -rf data/output data/output-demo data/fixture
