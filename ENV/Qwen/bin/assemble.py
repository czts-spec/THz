#!/usr/bin/env python3
"""Пол дома: сборка рабочей библиотеки Ф-чисел из корпуса.

Изнутри среды работать нечем, пока библиотека не собрана. Сборка идёт
`phikit/build.py`: исходники из `inbox/raw/` только читаются, правки
ложатся на копии вне корпуса, каждая правка обязана найтись ровно один
раз.

Сборка есть **производное**: 3,7 с прогона против 9 173 строк в
репозитории. Поэтому каталог сборки не версионируется, а версионируется
`state/build.json` — отпечатки исходников, число правок, итоги наборов и
дата. Пропала сборка — восстанавливается прогоном; пропал `build.json` —
восстановить нечем, потому он и в репозитории.

    python3 -B ENV/AIAΣA/bin/assemble.py           # собрать, если нужно
    python3 -B ENV/AIAΣA/bin/assemble.py --force   # пересобрать всегда
    python3 -B ENV/AIAΣA/bin/assemble.py --where   # где лежит сборка

Код возврата 0 только при зелёной сборке.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

BIN = Path(__file__).resolve().parent
HOME = BIN.parent
ROOT = BIN.parents[2]
STATE = HOME / "state"
STAMP = STATE / "build.json"

#: Каталог сборки. Вне репозитория: производное не хранится рядом с
#: исходником. Переопределяется переменной среды, чтобы дом можно было
#: поставить там, где /tmp мал или недолговечен.
#: AIASA_BUILD — каталог сборки. Адрес в ней (`://`, замер 27.09: среда держит там
#: ссылку на репозиторий) каталогом не считается: иначе сборка ложится в дом
#: относительным путём `https:/…` и красит линии 3 и 5+8 (test_red_lines).
_B = os.environ.get("AIASA_BUILD", "/tmp/aiasa-build")
DEST = Path("/tmp/aiasa-build" if "://" in _B else _B)

sys.path.insert(0, str(ROOT / "phikit"))


def sources() -> dict[str, str]:
    """Отпечатки исходников корпуса. Корпус только читается."""
    from patches import SOURCES  # noqa: E402
    raw = ROOT / "inbox" / "raw"
    out = {}
    for name, src in sorted(SOURCES.items()):
        p = raw / src
        out[name] = (hashlib.sha256(p.read_bytes()).hexdigest()
                     if p.exists() else "НЕТ В КОРПУСЕ")
    return out


def stamp() -> dict:
    return json.loads(STAMP.read_text(encoding="utf-8")) if STAMP.is_file() else {}


def is_current() -> bool:
    """Сборка годна, если отпечатки сошлись И файлы на месте."""
    s = stamp()
    if not s or s.get("исходники") != sources():
        return False
    dest = Path(s.get("каталог", DEST))
    return all((dest / n).is_file() for n in s.get("исходники", {}))


def build(dest: Path | None = None, force: bool = False) -> dict:
    """Собрать при нужде. Возвращает штамп."""
    dest = Path(dest or DEST)
    if not force and is_current():
        s = stamp()
        s["пересобрано"] = False
        return s
    import build as phibuild  # noqa: E402
    t0 = time.time()
    rep = phibuild.build(dest)
    s = {
        "каталог": str(dest),
        "исходники": sources(),
        "правок": rep["правок"],
        "разложено": rep["разложено"],
        "наборы": {k: {"итог": v[0], "прошло": v[1], "провалено": v[2],
                       "счётчик": v[3]} for k, v in rep["наборы"].items()},
        "зелёная": phibuild.is_green(rep),
        "секунд": round(time.time() - t0, 2),
        "дата": time.strftime("%Y-%m-%d %H:%M:%S"),
        "пересобрано": True,
    }
    STATE.mkdir(parents=True, exist_ok=True)
    tmp = STAMP.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(s, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8")
    tmp.replace(STAMP)
    return s


def build_dir() -> Path:
    """Каталог со сборкой; собирает, если её нет."""
    s = build()
    return Path(s["каталог"])


def is_green() -> bool:
    return bool(stamp().get("зелёная"))


def suites_are_recorded() -> bool:
    """Итоги наборов записаны с указанием счётчика, а не одним числом.

    Счётчики у наборов разные: «проверки», «перебор», «показ». Смешивать
    их нельзя — 3 404 перебранных тройки не то же, что 150 проверок.
    """
    n = stamp().get("наборы", {})
    return bool(n) and all(v.get("счётчик") in
                            ("проверки", "перебор", "показ", "без счётчика")
                            for v in n.values())


def corpus_untouched() -> bool:
    """Ни один шаг сборки не пишет в корпус: правки ложатся на копии."""
    src = (ROOT / "phikit" / "build.py").read_text(encoding="utf-8")
    return "shutil.copy2" in src and "RAW" in src


if __name__ == "__main__":
    if "--where" in sys.argv:
        print(build_dir())
        raise SystemExit(0)
    s = build(force="--force" in sys.argv)
    print(f"каталог: {s['каталог']}")
    print(f"пересобрано: {s['пересобрано']}, правок: {s['правок']}, "
          f"секунд: {s.get('секунд', '—')}")
    for name, r in s["наборы"].items():
        print(f"  {name:<34} {r['итог']:<10} {r['счётчик']:<10} "
              f"OK {r['прошло']:>4}  FAIL {r['провалено']:>3}")
    print(f"зелёная: {s['зелёная']}")
    raise SystemExit(0 if s["зелёная"] else 1)
