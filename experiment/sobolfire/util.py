"""Small helpers: checksums, atomic writes, timestamps."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


class ConfigError(ValueError):
    """A problem with inputs or settings; reported without a traceback."""


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_text(path, text):
    """Write via a temporary file and rename, so readers never see a partial file."""
    path = Path(path)
    tmp = path.with_name(f".{path.name}.partial")
    tmp.write_text(text)
    os.replace(tmp, path)


def write_json(path, obj):
    write_text(path, json.dumps(obj, indent=1, default=str) + "\n")


def read_json(path):
    with open(path) as f:
        return json.load(f)
