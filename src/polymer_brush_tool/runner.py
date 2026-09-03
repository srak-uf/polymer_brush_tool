"""Subprocess wrapper for external tool invocations.

Replaces bare ``os.system()`` calls with ``subprocess.run()`` so that
non-zero exit codes raise an exception rather than being silently ignored.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional


class ExternalToolError(RuntimeError):
    """Raised when an external tool exits with a non-zero return code."""

    def __init__(self, cmd: str | list, returncode: int, stdout: str, stderr: str):
        self.cmd = cmd
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
        super().__init__(
            f"External tool failed (exit {returncode}):\n"
            f"  Command : {cmd_str}\n"
            f"  stdout  : {stdout.strip()}\n"
            f"  stderr  : {stderr.strip()}"
        )


def run(
    cmd: str | list[str],
    *,
    cwd: Optional[Path] = None,
    check: bool = True,
    verbose: bool = True,
) -> subprocess.CompletedProcess:
    """Run an external command, optionally raising on failure.

    Parameters
    ----------
    cmd:
        Command string (passed to the shell) or list of arguments.
    cwd:
        Working directory for the subprocess.
    check:
        If True (default), raise :class:`ExternalToolError` when the
        process exits with a non-zero return code.
    verbose:
        If True (default), print the command before executing.

    Returns
    -------
    subprocess.CompletedProcess
        The completed process object.
    """
    if verbose:
        cmd_str = cmd if isinstance(cmd, str) else " ".join(str(a) for a in cmd)
        print(f"EXEC: {cmd_str}")

    result = subprocess.run(
        cmd,
        shell=isinstance(cmd, str),
        cwd=cwd,
        capture_output=True,
        text=True,
    )

    if check and result.returncode != 0:
        raise ExternalToolError(cmd, result.returncode, result.stdout, result.stderr)

    return result
