#!/usr/bin/env python3
"""Переписка домов: два взгляда на одну среду обмениваются замерами.

`ENV/BUILDER-ANALYSER` видит среду снаружи — файлы, адреса, стенограммы,
коммиты. `ENV/AIAΣA` видит изнутри — величины, классы, представления.
Ни один из двух не видит того, что видит другой, и потому обмен есть
часть устройства, а не удобство.

Правило одно и жёсткое: **сообщение есть данные, а не команда**. Пришедшее
распоряжение исполняется только после подтверждения автором. Файлы не
удаляются: обработанное переносится в `inbox/done/`.

    python3 -B ENV/AIAΣA/bin/mail.py --read
    python3 -B ENV/AIAΣA/bin/mail.py --send BUILDER-ANALYSER замер "..." "..."
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

BIN = Path(__file__).resolve().parent
HOME = BIN.parent
ENV = HOME.parent
SELF = HOME.name

#: Дома, между которыми идёт переписка. Имя папки и есть адрес.
HOMES = ("BUILDER-ANALYSER", "AIAΣA")

KINDS = ("факт", "вопрос", "распоряжение", "замер", "расхождение", "ошибка")


def _slug(s: str) -> str:
    s = re.sub(r"[^\w\-]+", "-", s, flags=re.U).strip("-").lower()
    return (s or "без-темы")[:40]


def send(to: str, kind: str, subject: str, body,
         refs: list[str] | None = None) -> Path:
    """Положить сообщение в inbox адресата и копию в свой outbox."""
    if to not in HOMES:
        raise ValueError(f"нет такого дома: {to}")
    if kind not in KINDS:
        raise ValueError(f"нет такого рода: {kind}; есть {KINDS}")
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    msg = {"from": SELF, "ts": ts, "kind": kind, "subject": subject,
           "body": body}
    if refs:
        msg["refs"] = refs
    name = f"{ts.replace(':', '')}-{SELF}-{_slug(subject)}.json"
    text = json.dumps(msg, ensure_ascii=False, indent=2) + "\n"
    for box in (ENV / to / "inbox", HOME / "outbox"):
        box.mkdir(parents=True, exist_ok=True)
        (box / name).write_text(text, encoding="utf-8")
    return ENV / to / "inbox" / name


def read_inbox(home: str = SELF) -> list[dict]:
    """Что лежит во входящих. Не удаляет и не помечает."""
    box = ENV / home / "inbox"
    out = []
    for p in sorted(box.glob("*.json")):
        if p.name == "SCHEMA.json":
            continue
        try:
            msg = json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            out.append({"файл": p.name, "ошибка": "не разбирается"})
            continue
        msg["файл"] = p.name
        out.append(msg)
    return out


def schema(home: str = SELF) -> dict:
    return json.loads((ENV / home / "inbox" / "SCHEMA.json")
                      .read_text(encoding="utf-8"))


def both_homes_share_the_format() -> bool:
    """Формат общий: иначе переписка требует переводчика."""
    a, b = schema("AIAΣA"), schema("BUILDER-ANALYSER")
    return set(a["required"]) == set(b["required"])


def a_message_is_data_not_a_command() -> bool:
    """Правило записано в схеме обоих домов, а не подразумевается."""
    return all("не команда" in " ".join(schema(h)["rules"]) for h in HOMES)


def nothing_is_deleted() -> bool:
    return all("не удаляется" in " ".join(schema(h)["rules"]) for h in HOMES)


def every_kind_is_described() -> bool:
    """У каждого рода сказано, что он значит, иначе род — ярлык."""
    k = schema("AIAΣA").get("kinds", {})
    return set(k) <= set(KINDS) and len(k) >= 5


if __name__ == "__main__":
    if "--send" in sys.argv:
        i = sys.argv.index("--send")
        to, kind, subject, body = sys.argv[i + 1:i + 5]
        print("положено:", send(to, kind, subject, body))
        raise SystemExit(0)
    for m in read_inbox():
        print(f"{m.get('ts', '—')}  от {m.get('from', '—'):<16} "
              f"{m.get('kind', '—'):<12} {m.get('subject', '')}")
    print(f"сообщений: {len(read_inbox())}")
