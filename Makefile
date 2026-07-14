.PHONY: setup data evaluate test clean all

all: data evaluate

setup:
	pip install -r requirements.txt

data:
	python src/panel_eval/simulate_counts.py

evaluate:
	python scripts/run_evaluation.py

test:
	pytest -q

clean:
	rm -rf results/figures/* results/tables/* results/REPORT.md
