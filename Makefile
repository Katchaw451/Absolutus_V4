.PHONY: install test bench
install:
	python3 -m venv .venv && . .venv/bin/activate && pip install -e ".[dev]"
test:
	python3 -m pytest -q
bench:
	python3 -m absolutus_v4.bench --seeds 20 --out results --docs docs/RESULTS_SIM.md
