"""Shared set-up for the load-model tests. The repository root is put on the
import path by the pytest setting in pyproject.toml.
Author: Akosa Samuel Onyejekwe (independent)"""
import os

import pytest


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """With UNISTALL_NO_SKIPS set (the continuous-integration run), a skipped
    test is reported as failed: there a skip means a check did not run."""
    outcome = yield
    report = outcome.get_result()
    if os.environ.get("UNISTALL_NO_SKIPS") and report.skipped:
        report.outcome = "failed"
        report.longrepr = f"skipped where no skip is allowed: {report.longrepr}"
