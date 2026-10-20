#!/usr/bin/env python3
"""autofix: выскочившее в чат ложится в реестр — устройством, не дисциплиной.

Слово автора 02.09: «переписка в чате — дурной тон. в чат только ссылки…
ответ — это вопрос протокола и дисциплины… индексированный реестр —
это скрипт, против вспоминания токенов. дополни протокол обмена так,
чтобы выскочившее в чат так же ложилось в реестр, с пометкой autofix…
это п��рерасход дорогих токенов, то есть ситуация аварийная… ставь
крючок на автофикс, и в процедуре фикса напоминай себе протокол».

Что делает: читает живую стенограмму сессии (профиль Claude Code —
тот же источник, что у tools/export_chat.py), находит реплики
исполнителя, ушедшие в чат ВОПРЕКИ протоколу (не ссылки, а изложение),
и кладёт их в реестр разделами с пометкой autofix — по дням, с
ограничением размера части. Печатает напоминание протокола. Водяной
знак — по времени последней обработанной реплики, на сессию.

    python3 -B ENV/AIAΣA/bin/autofix.py          # положить и отчитаться
    python3 -B ENV/AIAΣA/bin/autofix.py --dry    # только показать
    python3 -B ENV/AIAΣA/bin/autofix.py --hook   # из крючка: коротко + напоминание

Ядро принимает пути параметрами — приёмка гоняет его на песочнице.
"""
from __future__ import annotations

import json
import pathlib
import re
import sys
import urllib.parse

NOT_MEASURED = (
    "порог изложения (200 знаков без ссылок) — граница, не замер: короткая "
    "реплика по существу пройдёт мимо, длинный перечень ссылок с подписями — "
    "тоже; что именно потеряно, не считается",
    "стенограмма читается из профиля Claude Code этой машины; чужие сессии "
    "и чат вне профиля инструменту не видны",
)

BIN = pathlib.Path(__file__).resolve().parent
HOME = BIN.parent
ROOT = HOME.parent.parent
PROJECTS = pathlib.Path.home() / ".claude" / "projects"
CHATLOG = ROOT / "chatlog" / "raw"
REGISTER = ROOT / "notes" / "регистр"
INDEX = ROOT / "notes" / "РЕГИСТР внешней сверки.md"

MIN_PROSE = 200
PART_CAP = 6_500

PROTOCOL = ("ПРОТОКОЛ: в чат — только ссылки, собранные машинно "
            "(python3 -B ENV/AIAΣA/bin/link.py <путь>); ответ — файлом в "
            "notes/ и коммит; пункты автора — обязательная очередь, не меню. "
            "Ушедшее в чат — перерасход дорогих токенов, аварийный режим; "
            "autofix кладёт его в реестр, но не отменяет правила.")

_STRIP = (re.compile(r"`[^`]*`"), re.compile(r"https?://\S+"),
          re.compile(r"\S*[/\\]\S*"), re.compile(r"\b[0-9a-f]{7,40}\b"))
_HEADING = re.compile(r"^## §?(\d+)([а-я])?\.?\s")
LINK_TAIL = 40


def is_link_line(s: str) -> bool:
    rest = s
    for rx in _STRIP:
        rest = rx.sub("", rest)
    rest = rest.strip(" -*•:;,.—()")
    return len(rest) <= LINK_TAIL and rest != s.strip()


def prose_length(text: str) -> int:
    n, in_code = 0, False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("```"):
            in_code = not in_code
            continue
        if in_code or not s or s.startswith("#") or is_link_line(s):
            continue
        n += len(s)
    return n


def messages(files) -> list[dict]:
    out = []
    for fn in files:
        sid = re.sub(r"\.part-\d+$", "", pathlib.Path(fn).stem)
        try:
            fh = open(fn, encoding="utf-8", errors="ignore")
        except OSError:
            continue
        with fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue
                if rec.get("type") != "assistant":
                    continue
                content = rec.get("message", {}).get("content")
                if not isinstance(content, list):
                    continue
                text = "\n\n".join(c.get("text", "") for c in content
                                    if isinstance(c, dict)
                                    and c.get("type") == "text"
                                    and c.get("text", "").strip())
                if not text:
                    continue
                out.append({"session": sid, "ts": rec.get("timestamp", ""),
                            "uuid": rec.get("uuid", ""), "text": text})
    uniq = {}
    for m in out:
        uniq.setdefault(m["uuid"] or (m["ts"], m["text"][:40]), m)
    out = list(uniq.values())
    out.sort(key=lambda m: (m["ts"], m["uuid"]))
    return out


_ITEM = re.compile(r"^### \S+ · (\S+)\s*$")


def recorded(register: pathlib.Path) -> set[str]:
    seen = set()
    for p in register.glob("*autofix*.md"):
        for line in p.read_text(encoding="utf-8").splitlines():
            m = _ITEM.match(line)
            if m:
                seen.add(m.group(1))
    return seen


def pending(msgs, seen: set[str]) -> list[dict]:
    return [m for m in msgs
            if m["uuid"] and m["uuid"] not in seen
            and prose_length(m["text"]) >= MIN_PROSE]


def _remote_labels(register: pathlib.Path) -> list[int]:
    import subprocess
    root = register.parents[1]
    r = subprocess.run(["git", "ls-tree", "--name-only", "origin/main",
                        "notes/регистр/"],
                       cwd=root, capture_output=True, text=True)
    if r.returncode != 0:
        return []
    nums = []
    for line in r.stdout.splitlines():
        m = re.search(r"/(\d{3}) ", line)
        if m:
            nums.append(int(m.group(1)))
    return nums


def _labels(register: pathlib.Path) -> list[int]:
    nums = []
    for p in register.glob("*.md"):
        m = re.match(r"^(\d{3})", p.name)
        if m:
            nums.append(int(m.group(1)))
    return sorted(nums + _remote_labels(register))


def _sanitize(text: str) -> str:
    return "\n".join(("▸ " + l) if l.startswith("## ") else l
                      for l in text.splitlines())


def _fmt_size(n: int) -> str:
    return f"{n:,}".replace(",", " ")


def _index_row(label: int, title: str, name: str, size: int) -> str:
    link = "регистр/" + urllib.parse.quote(name)
    cell = title.replace("|", "/").replace(",", "")
    return f"| §{label} | [{cell}]({link}) | {_fmt_size(size)} |"


def _write_section(register, index, label, title, body) -> pathlib.Path:
    name = f"{label:03d} — {title}.md"
    path = register / name
    text = f"## §{label}. {title}\n\n{body.rstrip()}\n"
    path.write_text(text, encoding="utf-8")
    row = _index_row(label, title, name, len(text.encode("utf-8")))
    idx = index.read_text(encoding="utf-8")
    if not idx.endswith("\n"):
        idx += "\n"
    index.write_text(idx + row + "\n", encoding="utf-8")
    return path


def index_fill(register=None, index=None) -> dict:
    register = register or REGISTER
    index = index or INDEX
    текст = index.read_text(encoding="utf-8")
    строки = текст.splitlines()
    добавлено = []
    снято = []
    живые = []
    for l in строки:
        m_link = re.search(r"\]\(регистр/([^)]+)\)", l)
        if m_link and not (register / urllib.parse.unquote(m_link.group(1))).exists():
            снято.append(l)
            continue
        живые.append(l)
    if снято:
        строки = живые
        текст = "\n".join(строки) + "\n"
    for part in sorted(register.glob("*.md")):
        m = re.match(r"^(\d+) — (.+)\.md$", part.name)
        if not m:
            continue
        link = "регистр/" + urllib.parse.quote(part.name)
        if link in текст:
            continue
        заголовок = part.read_text(encoding="utf-8").splitlines()[0]
        загол = re.sub(r"^##\s*§\d+\.\s*", "", заголовок).strip()
        row = _index_row(int(m.group(1)), загол or m.group(2), part.name,
                         part.stat().st_size)
        место = len(строки)
        for i, l in enumerate(строки):
            n = re.match(r"\| §(\d+) \|", l)
            if n and int(n.group(1)) > int(m.group(1)):
                место = i
                break
        строки.insert(место, row)
        добавлено.append(part.name)
    if добавлено or снято:
        index.write_text("\n".join(строки) + "\n", encoding="utf-8")
    return {"досыпано строк": len(добавлено), "части": добавлено, "снято строк": len(снято)}


_COUNT = re.compile(r"реплик: (\d+)\.")


def _day_section(register, day):
    cands = sorted(register.glob(f"*autofix: в чат вопреки протоколу, {day}*.md"))
    return cands[-1] if cands else None


def _append_items(register, index, path, items) -> None:
    text = path.read_text(encoding="utf-8")
    m = _COUNT.search(text)
    n = int(m.group(1)) if m else 0
    text = _COUNT.sub(f"реплик: {n + len(items)}.", text, count=1)
    text = text.rstrip() + "\n\n" + "".join(items)
    path.write_text(text, encoding="utf-8")
    size = _fmt_size(len(text.encode("utf-8")))
    link = "регистр/" + urllib.parse.quote(path.name)
    rows = index.read_text(encoding="utf-8").split("\n")
    rows = [re.sub(r"\| [\d ]+ \|\s*$", f"| {size} |", r) if link in r else r for r in rows]
    index.write_text("\n".join(rows), encoding="utf-8")


def record(entries, *, register=REGISTER, index=INDEX) -> list[pathlib.Path]:
    if not entries:
        return []
    by_day: dict[str, list] = {}
    for m in entries:
        by_day.setdefault(m["ts"][:10] or "без-даты", []).append(m)
    written = []
    nxt = (max(_labels(register)) if _labels(register) else 0) + 1
    for day in sorted(by_day):
        rest = list(by_day[day])
        existing = _day_section(register, day)
        if existing is not None:
            have = len(existing.read_text(encoding="utf-8"))
            take = []
            while rest:
                item = (f"### {rest[0]['ts']} · {rest[0]['uuid']}\n\n"
                        f"{_sanitize(rest[0]['text'])}\n\n---\n")
                if have + len(item) > PART_CAP:
                    break
                take.append(item); have += len(item); rest.pop(0)
            if take:
                _append_items(register, index, existing, take)
                written.append(existing)
        if not rest:
            continue
        chunks, cur, size = [], [], 0
        for m in rest:
            item = (f"### {m['ts']} · {m['uuid']}\n\n"
                    f"{_sanitize(m['text'])}\n\n---\n")
            if cur and size + len(item) > PART_CAP:
                chunks.append(cur)
                cur, size = [], 0
            cur.append(item)
            size += len(item)
        if cur:
            chunks.append(cur)
        for i, chunk in enumerate(chunks, 1):
            suffix = f" ({i} из {len(chunks)})" if len(chunks) > 1 else ""
            title = f"autofix: в чат вопреки протоколу, {day}{suffix}"
            head = (f"`[вычислено здесь]` **autofix.** Реплики исполнителя, "
                    f"ушедшие в чат вопреки протоколу (в чат — только "
                    f"ссылки; ответ — файлом). Аварийный режим: перерасход "
                    f"дорогих токенов. Положено в реестр устройством "
                    f"(`ENV/AIAΣA/bin/autofix.py`), дословно, сессия "
                    f"`{by_day[day][0]['session']}`; реплик: {len(chunk)}.\n\n"
                    f"{PROTOCOL}\n\n")
            written.append(_write_section(register, index, nxt, title,
                                          head + "".join(chunk)))
            nxt += 1
    return written


def run(*, files=None, register=REGISTER, index=INDEX, dry=False):
    if files is None:
        files = sorted(set(list(PROJECTS.glob("*/*.jsonl"))
                            + list(PROJECTS.glob("*.jsonl"))
                            + list(CHATLOG.glob("*.part-*.jsonl"))))
    msgs = messages(files)
    todo = pending(msgs, recorded(register))
    written = [] if dry else record(todo, register=register, index=index)
    fill = {"досыпано строк": 0, "части": []} if dry else index_fill(register, index)
    return {"реплик": len(msgs), "ушло_в_чат": len(todo),
            "разделов": [p.name for p in written],
            "досыпано строк оглавления": fill["досыпано строк"]}


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    try:
        r = run(dry=(mode == "--dry"))
    except Exception as ex:
        print(f"autofix: сбой {type(ex).__name__}: {ex}")
        sys.exit(0)
    if mode == "--hook":
        if r["разделов"]:
            print(f"autofix: {r['ушло_в_чат']} реплик легли в реестр: "
                  f"{', '.join(r['разделов'])}")
        print(PROTOCOL)
    elif mode == "--dry":
        print(f"реплик исполнителя: {r['реплик']}; ушло в чат (не записано): "
              f"{r['ушло_в_чат']}")
    else:
        print(f"реплик: {r['реплик']}; ушло в чат: {r['ушло_в_чат']}; "
              f"разделов положено: {len(r['разделов'])}; "
              f"строк оглавления досыпано: {r['досыпано строк оглавления']}")
        for n in r["разделов"]:
            print("  ", n)
    sys.exit(0)
