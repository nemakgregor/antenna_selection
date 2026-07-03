import json
import os


def atomic_write_csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".tmp")
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(tmp_path, index=False, compression=compression)
    os.replace(tmp_path, path)


def atomic_write_json(payload, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".tmp")
    tmp_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(tmp_path, path)
