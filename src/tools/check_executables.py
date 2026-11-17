"""Check the affinity executables required by validated pipeline options."""

import json
from pathlib import Path
import shutil
import subprocess
import sys


def check_executables(executables):
    for name in executables:
        path = shutil.which(name)
        if not path:
            raise ValueError(f"required executable '{name}' was not found in PATH")
        try:
            probe = subprocess.run([path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, errors="replace", timeout=15, check=False)
        except (OSError, subprocess.TimeoutExpired) as error:
            raise ValueError(f"required executable '{name}' could not be checked: {error}") from error
        output = probe.stdout or ""
        failures = ("no binaries found", "bad cpu type", "cannot execute", "exec format error",
                    "error while loading shared libraries")
        if probe.returncode < 0 or probe.returncode in (126, 127, 132, 134, 139) or any(
            message in output.lower() for message in failures
        ):
            raise ValueError(f"required executable '{name}' is present but cannot run:\n{output.strip()}")


if __name__ == "__main__":
    try:
        options = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
        check_executables(options["required_executables"])
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"ERROR: executable preflight failed: {error}")
