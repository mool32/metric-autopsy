"""Shared pytest configuration.

Every run with a pre-registration writes to the run log by default. Point each test at its
own log file so the suite never writes into the working tree and attempt counters start at 1.
"""
import pytest


@pytest.fixture(autouse=True)
def _isolated_run_log(tmp_path, monkeypatch):
    monkeypatch.setenv("METRIC_AUTOPSY_LOG", str(tmp_path / "metric_autopsy_runs.jsonl"))
