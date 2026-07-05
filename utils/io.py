import json
import os
import tarfile


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


def archive_csv_artifacts(directory, archive_name="csv_data.tar.gz", remove_originals=True):
    directory = os.fspath(directory)
    archive_path = os.path.join(directory, archive_name)
    csv_paths = []
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if not os.path.isfile(path):
            continue
        if name == archive_name:
            continue
        if name.endswith(".csv") or name.endswith(".csv.gz"):
            csv_paths.append(path)

    if not csv_paths:
        return archive_path

    tmp_path = archive_path + ".tmp"
    with tarfile.open(tmp_path, "w:gz") as archive:
        for path in csv_paths:
            archive.add(path, arcname=os.path.basename(path))
    os.replace(tmp_path, archive_path)

    if remove_originals:
        for path in csv_paths:
            os.remove(path)

    return archive_path
