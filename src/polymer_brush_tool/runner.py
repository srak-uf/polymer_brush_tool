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

    def __init__(self, cmd: str | list, returncode: int, stdout: str = "", stderr: str = ""):
        self.cmd = cmd
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        cmd_str = cmd if isinstance(cmd, str) else " ".join(str(a) for a in cmd)
        msg = f"External tool failed (exit {returncode}):\n  Command : {cmd_str}"
        if stdout.strip():
            msg += f"\n  stdout  : {stdout.strip()[-2000:]}"
        if stderr.strip():
            msg += f"\n  stderr  : {stderr.strip()[-2000:]}"
        super().__init__(msg)


def run(
    cmd: str | list[str],
    *,
    cwd: Optional[Path] = None,
    check: bool = True,
    verbose: bool = True,
    capture: bool = False,
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
    capture:
        If True, capture stdout/stderr and attach them to the returned
        object (and to the error).  Default False so that long-running
        tools such as ``gmx mdrun -v`` stream their progress to the
        terminal exactly as they did with ``os.system``.

    Returns
    -------
    subprocess.CompletedProcess
    """
    if verbose:
        cmd_str = cmd if isinstance(cmd, str) else " ".join(str(a) for a in cmd)
        print(f"EXEC: {cmd_str}", flush=True)

    result = subprocess.run(
        cmd,
        shell=isinstance(cmd, str),
        cwd=cwd,
        capture_output=capture,
        text=True,
    )

    if check and result.returncode != 0:
        raise ExternalToolError(
            cmd, result.returncode, result.stdout or "", result.stderr or ""
        )

    return result
