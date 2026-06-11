import numpy as np
import pandas as pd


def format_number(value, precision=4):
    if pd.isna(value):
        return ""
    if abs(value) >= 1e5 or (0 < abs(value) < 1e-3):
        return f"{value:.{precision}e}"
    return f"{value:.{precision}f}"


def format_sigma(value):
    return f"{value:g}"


def format_number_slug(value):
    if float(value).is_integer():
        return str(int(value))
    return str(value).replace(".", "p")


def write_markdown(path, lines):
    path.write_text("\n".join(lines), encoding="utf-8")


def leader_segments(rows, x_col="sigma", leader_col="leader"):
    segments = []
    current = None
    start_x = None
    end_x = None
    for _, row in rows.sort_values(x_col).iterrows():
        leader = row[leader_col]
        x_value = row[x_col]
        if current is None:
            current = leader
            start_x = x_value
            end_x = x_value
            continue
        if leader == current:
            end_x = x_value
            continue
        segments.append((start_x, end_x, current))
        current = leader
        start_x = x_value
        end_x = x_value
    if current is not None:
        segments.append((start_x, end_x, current))
    return segments


def format_leader_segments(segments, format_x=str):
    return ", ".join(
        (
            f"{format_x(start)}: {leader}"
            if start == end
            else f"{format_x(start)}..{format_x(end)}: {leader}"
        )
        for start, end, leader in segments
    )


def split_win_shares(
    runs,
    case_cols,
    winner_group_cols,
    sample_group_cols,
    algorithm_col,
    metric_resolver,
    sample_col="sample",
    sort_cols=None,
):
    rows = []
    for keys, chunk in runs.groupby(case_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        case = dict(zip(case_cols, keys))
        by_algorithm = chunk.set_index(algorithm_col)
        for metric, label, direction in metric_resolver(case, chunk):
            values = by_algorithm[metric]
            best_value = values.min() if direction == "min" else values.max()
            winners = values[np.isclose(values, best_value)].index.tolist()
            share = 1.0 / len(winners)
            for winner in winners:
                rows.append(
                    {
                        **case,
                        "metric": metric,
                        "metric_label": label,
                        algorithm_col: winner,
                        "win_share": share,
                        "winner_hit": 1.0,
                    }
                )

    wins = pd.DataFrame(rows)
    group_cols = [*winner_group_cols, "metric", "metric_label", algorithm_col]
    result = (
        wins.groupby(group_cols, as_index=False)
        .agg(win_share=("win_share", "sum"), winner_hits=("winner_hit", "sum"))
        .merge(
            runs.groupby(sample_group_cols, as_index=False)[sample_col]
            .nunique()
            .rename(columns={sample_col: "samples"}),
            on=sample_group_cols,
            how="left",
        )
        .assign(
            win_fraction=lambda df: df["win_share"] / df["samples"],
            winner_rate=lambda df: df["winner_hits"] / df["samples"],
        )
    )
    return result.sort_values(sort_cols) if sort_cols else result
