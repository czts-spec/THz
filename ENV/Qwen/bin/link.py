#!/usr/bin/env python3
"""Ссылка на файл репозитория — собирается машинно, не пишется рукой.

Слово автора 03.09: «.md должны быть ссылками, простое упоминание
неприемлемо — я это всё читаю». Правило §241 (register_split.url):
адрес не сочиняется, а считается из имени файла на диске; вручную
набранный адрес однажды дал 404.

    python3 -B ENV/AIAΣA/bin/link.py <путь…>        # markdown-ссылки
    python3 -B ENV/AIAΣA/bin/link.py --url <путь>   # голый адрес
"""
from __future__ import annotations

import pathlib
import subprocess
import sys
import urllib.parse

ROOT = pathlib.Path(__file__).resolve().parents[3]
BLOB = "https://github.com/dtba3a-del/InvesePolar/blob/main/"
#: Ветка, по которой открывается ссылка: файл обязан лежать В НЕЙ, не на диске.
REMOTE = "origin/main"

#: Слово автора 17.09 (коррекция `bb2912c048`): ссылка на ответ ушла в чат за
#: 285 минут до пуша — файл лежал только в рабочем дереве, и автор получил 404.
#: Прежняя пометка NOT_MEASURED этого места («доставка в origin — нет») была
#: снята замером в ту же ночь: устройство спрашивает у git, лежит ли путь в
#: origin/main, и без этого ссылку НЕ собирает. `--local` — осознанный обход
#: для служебного текста, который сам говорит, что файл ещё едет.
NOT_MEASURED = ()


class NotDelivered(FileNotFoundError):
    """Файл есть на диске, но в origin/main его нет: ссылка дала бы 404."""


def delivered(rel: str, remote: str = REMOTE) -> bool:
    """`[вычислено здесь]` Путь лежит в удалённой ветке (git cat-file -e)."""
    r = subprocess.run(["git", "-C", str(ROOT), "cat-file", "-e", f"{remote}:{rel}"],
                       capture_output=True)
    return r.returncode == 0


def url(path, local: bool = False) -> str:
    p = pathlib.Path(path)
    rel = (p.resolve().relative_to(ROOT) if p.is_absolute()
           else p).as_posix()
    if not (ROOT / rel).exists():
        raise FileNotFoundError(f"нет такого файла: {rel}")
    if not local and not delivered(rel):
        raise NotDelivered(f"файл не запушен — в {REMOTE} его нет, ссылка дала бы 404: {rel}")
    return BLOB + urllib.parse.quote(rel)


def md(path, label: str | None = None, local: bool = False) -> str:
    rel = pathlib.Path(path).as_posix()
    name = label or pathlib.Path(rel).name
    return f"[{name}]({url(rel, local=local)})"


if __name__ == "__main__":
    args = sys.argv[1:]
    local = "--local" in args
    args = [a for a in args if a != "--local"]
    bare = args and args[0] == "--url"
    for a in (args[1:] if bare else args):
        try:
            print(url(a, local=local) if bare else md(a, local=local))
        except FileNotFoundError as ex:
            print(f"ОТКАЗ: {ex}")
            sys.exit(1)
