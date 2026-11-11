"""independence.py: a session's transcript against the paths it may not name (no network, no real transcript)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import independence  # noqa: E402


def _line(role, block):
    return json.dumps({"type": role, "message": {"role": role, "content": [block]}})


def _transcript(calls, results):
    out = []
    for c, r in zip(calls, results):
        out.append(_line("assistant", {"type": "tool_use", "name": "Read", "input": {"file_path": c}}))
        out.append(_line("user", {"type": "tool_result", "content": r}))
    return "\n".join(out) + "\n"


def _setup(tmp_path, transcripts):
    tr = tmp_path / "tasks"
    tr.mkdir()
    agents = []
    for i, (p, b, text) in enumerate(transcripts):
        sid = f"s{i}"
        (tr / f"{sid}.output").write_text(text)
        agents.append(f"{p}\t{b}\t{sid}")
    (tmp_path / "agents.tsv").write_text("\n".join(agents) + "\n")
    for p in ("1", "2"):
        d = tmp_path / f"pass{p}"
        d.mkdir()
        (d / "batch_01.list").write_text("R00001\ttexts/R00001.txt\n")
    return tr


def test_calls_only_counts_the_session_s_own_requests(tmp_path, capsys):
    instructions = "the protocol (validation/plans/level_b_protocol.md) is cited here"
    reads_it = _transcript([f"{tmp_path}/pass1/batch_01.list", "validation/plans/level_b_protocol.md"],
                           ["R00001\ttexts/R00001.txt", "the protocol's text"])
    cites_it = _transcript([f"{tmp_path}/pass2/batch_01.list", "validation/level_b/screening.md"],
                           ["R00001\ttexts/R00001.txt", instructions])
    tr = _setup(tmp_path, [("1", "01", reads_it), ("2", "01", cites_it)])
    args = ["--transcripts", str(tr), "--agents", str(tmp_path / "agents.tsv"),
            "--pass1", f"{tmp_path}/pass1", "--pass2", f"{tmp_path}/pass2",
            "--lists", f"{tmp_path}/pass1", f"{tmp_path}/pass2", "--batch-pattern", r"/batch_(\d+)\.list"]
    # as a forbidden path, the cited protocol fails both sessions
    assert independence.main(args + ["--forbid", "level_b_protocol"]) == 1
    out = capsys.readouterr().out
    assert out.count("CHECK") == 2
    # as a path only the session's own calls may not name, only the session that opened it fails
    assert independence.main(args + ["--calls-only", "level_b_protocol"]) == 1
    lines = capsys.readouterr().out.splitlines()
    assert "in its own calls 1" in lines[0] and lines[0].endswith("CHECK")
    assert "level_b_protocol named 1 times, in its own calls 0" in lines[1] and lines[1].endswith("ok")
    assert lines[2] == "2 sessions: 1 to check"


def test_calls_text_reads_tool_calls_only():
    s = _transcript(["a/b.md"], ["c/d.md"]) + "not json\n"
    assert independence.calls_text(s) == json.dumps({"file_path": "a/b.md"})


def test_s4_s_texts_by_their_keys(tmp_path, capsys):
    """S4's batch lists hold its texts' keys (`V<nnnnn>`) first; a session naming another batch's text fails."""
    own = _transcript([f"{tmp_path}/pass1/batch_01.list", "texts/V00001.txt"], ["V00001\tR00189:1\tMED:1\ttexts/V00001.txt", "a text"])
    other = _transcript([f"{tmp_path}/pass2/batch_01.list", "texts/V00002.txt"], ["V00001\tR00189:1\tMED:1\ttexts/V00001.txt", "a text"])
    tr = _setup(tmp_path, [("1", "01", own), ("2", "01", other)])
    for p in ("1", "2"):
        (tmp_path / f"pass{p}" / "batch_01.list").write_text("V00001\tR00189:1\tMED:1\ttexts/V00001.txt\n")
    args = ["--transcripts", str(tr), "--agents", str(tmp_path / "agents.tsv"),
            "--pass1", f"{tmp_path}/pass1", "--pass2", f"{tmp_path}/pass2",
            "--lists", f"{tmp_path}/pass1", f"{tmp_path}/pass2", "--batch-pattern", r"/batch_(\d+)\.list",
            "--text-pattern", r"\b(V\d{5})\.txt"]
    assert independence.main(args) == 1
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].endswith("ok") and "outside its list ['V00002']" in lines[1] and lines[1].endswith("CHECK")


def test_the_harness_s_attachment_is_counted_apart(tmp_path, capsys):
    """S4 (J31): the harness's task list in a session's context names a forbidden path; the session's own reading
    does not. Strict, it fails; with the attachment counted apart, only a session that names the path itself fails."""
    att = json.dumps({"type": "attachment", "attachment": {"type": "task_reminder", "content": [
        {"id": "1", "subject": "S3: decision sheet from screening/to_owner.tsv and claims.md"}]}})
    clean = att + "\n" + _transcript([f"{tmp_path}/pass1/batch_01.list"], ["R00001\ttexts/R00001.txt"])
    reads = att + "\n" + _transcript([f"{tmp_path}/pass2/batch_01.list", "validation/level_b/claims.md"],
                                     ["R00001\ttexts/R00001.txt", "the claims"])
    tr = _setup(tmp_path, [("1", "01", clean), ("2", "01", reads)])
    args = ["--transcripts", str(tr), "--agents", str(tmp_path / "agents.tsv"),
            "--pass1", f"{tmp_path}/pass1", "--pass2", f"{tmp_path}/pass2", "--forbid", "claims.md",
            "--batch-pattern", r"/batch_(\d+)\.list"]
    assert independence.main(args) == 1
    assert capsys.readouterr().out.count("CHECK") == 2
    out_file = tmp_path / "excluded.jsonl"
    assert independence.main(args + ["--exclude-attachment", "task_reminder", "--excluded-out", str(out_file)]) == 1
    lines = capsys.readouterr().out.splitlines()
    assert "claims.md named 0 times" in lines[0] and "excluded attachments 1 (1 distinct), naming: claims.md 1" in lines[0]
    assert lines[0].endswith("ok") and "claims.md named 1 times" in lines[1] and lines[1].endswith("CHECK")
    assert lines[2] == "2 sessions: 1 to check" and lines[3].endswith("1, written to excluded.jsonl")
    assert len(out_file.read_text().splitlines()) == 1 and "to_owner.tsv" in out_file.read_text()
