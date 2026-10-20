#!/usr/bin/env python3
"""Указатель писем автора (ENV/AIAΣA/notes/Q) и ответов на них — ребро для
каждого письма в графе ссылок.

Слово автора 04.09 (письмо 22:40): острова и сироты — не мусор, а
законченные ходы мысли; мера islands.py показала, что письма в Q/ —
острова, потому что ответы ссылаются на них датой и содержанием, а не
именем файла. Этот указатель собирается прогоном (производное, руками не
правится) и пишется в notes/: каждое письмо — строка с размером, датой,
найденными ответами (файлы notes/ и ENV/AIAΣA/notes/A, где встречается имя
письма, его дата вместе с ключевым словом заголовка, или хеш загрузки) и
dолей чтения из state/letters_read.json (доля ставится только после
чтения на 100 % с названным числом байт). [вычислено здесь]
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
Q = ROOT / "ENV/AIAΣA/notes/Q"
A = ROOT / "ENV/AIAΣA/notes/A"
NOTES = ROOT / "notes"
STATE = ROOT / "ENV/AIAΣA/state/letters_read.json"
OUT = NOTES / "письма автора — указатель: ответы и доля чтения (собирается прогоном).md"


def read_state() -> dict:
    try:
        return json.loads(STATE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def keywords(name: str) -> list[str]:
    body = re.sub(r"^\d{4}-\d\d-\d\d\s*—\s*", "", name)
    body = re.sub(r"\.\w+$", "", body)
    words = [w for w in re.split(r"[^\w]+", body) if len(w) >= 6]
    return words[:4]


def answers_for(letter: pathlib.Path, texts: dict) -> list[str]:
    name = letter.name
    date = name[:10]
    kws = keywords(name)
    hits = []
    for p, t in texts.items():
        if name in t or letter.stem in t:
            hits.append(p); continue
        if date in t and kws and sum(1 for k in kws if k in t) >= max(2, len(kws) - 1):
            hits.append(p)
    return sorted(str(p.relative_to(ROOT)) for p in hits)


def build() -> str:
    letters = sorted(p for p in Q.glob("*") if p.is_file() and p.name != "README.md")
    texts = {}
    for d in (NOTES, A):
        for p in d.glob("*.md"):
            if p == OUT:
                continue
            try:
                texts[p] = p.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                pass
    st = read_state()
    rows, no_answer, unread = [], 0, 0
    for p in letters:
        ans = answers_for(p, texts)
        rd = st.get(p.name, {})
        share = rd.get("доля", "—")
        if not ans:
            no_answer += 1
        if share == "—":
            unread += 1
        a_txt = "; ".join(f"`{a}`" for a in ans[:3]) + (" …" if len(ans) > 3 else "") if ans else "**ответа по имени/дате не найдено**"
        size = f"{p.stat().st_size:,}".replace(",", " ")
        rows.append(f"| `{p.name}` | {size} | {share} | {a_txt} |")
    head = [
        "# Письма автора — указатель: ответы и доля чтения (собирается прогоном)",
        "",
        "Производное: `python3 -B ENV/AIAΣA/bin/letters_index.py --write`; руками не",
        "правится. Основание — письмо автора 04.09 22:40 (острова — сокровище) и мера",
        "`islands.py`: письма без ребра к ответам были островами графа ссылок. Доля",
        "чтения — из `ENV/AIAΣA/state/letters_read.json`, ставится только после чтения",
        "на 100 % (байты названы). [вычислено здесь]",
        "",
        f"Писем: {len(letters)}; без найденного ответа: {no_answer}; без замера ч��ения: {unread}.",
        "",
        "| письмо | байт | прочитано | ответы (файлы, где найдено имя/дата+ключевые слова) |",
        "|---|---|---|---|",
    ]
    return "\n".join(head + rows) + "\n"


def main() -> int:
    text = build()
    if "--write" in sys.argv:
        OUT.write_text(text, encoding="utf-8")
        print(f"записано: {OUT.name}")
    print("\n".join(text.splitlines()[8:9]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
