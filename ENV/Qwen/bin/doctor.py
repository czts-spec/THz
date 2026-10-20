#!/usr/bin/env python3
"""Самопроверка дома AIAΣA. Точка входа после перерыва.

Код возврата 0 — дом цел, библиотека собрана, инструменты отвечают. Сети
не требует. Первый запуск дольше остальных: он собирает библиотеку.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HOME = Path(__file__).resolve().parent.parent
BIN = HOME / "bin"
ROOT = HOME.parents[1]
sys.path.insert(0, str(BIN))

REQUIRED_DIRS = ("bin", "state", "inbox", "outbox", "notes", "skills")
REQUIRED_FILES = ("MANIFEST.json", "RUNBOOK.json", "env.json",
                   "inbox/SCHEMA.json", "state/build.json", "CLAUDE.md")
TOOLS = ("assemble.py", "phi.py", "classes.py", "arith.py", "translate.py",
          "dual.py", "probe_env.py", "mail.py", "doctor.py")
SKILLS = ("phi-first", "lazy-translation", "dual-path", "named-denominator",
           "gap-register", "false-class-guard")
NOTES = ("00-ПОРЯДОК.md", "10-КЛАССЫ.md", "20-ОБМЕН.md", "30-ПУСТОТЫ.md",
          "40-ИСТОЧНИКИ.md", "50-ВХОД.md")

FAIL: list[str] = []
OK: list[str] = []


def check(name: str, cond: bool) -> None:
    (OK if cond else FAIL).append(name)


def main() -> int:
    for d in REQUIRED_DIRS:
        check(f"каталог {d}", (HOME / d).is_dir())
    for f in REQUIRED_FILES:
        check(f"файл {f}", (HOME / f).is_file())
    for t in TOOLS:
        check(f"инструмент {t}", (BIN / t).is_file())
    for s in SKILLS:
        check(f"скил {s}", (HOME / "skills" / s / "SKILL.md").is_file())
    for n in NOTES:
        check(f"памятка {n}", (HOME / "notes" / n).is_file())

    for f in ("MANIFEST.json", "RUNBOOK.json", "env.json",
              "inbox/SCHEMA.json", "state/build.json"):  # CLAUDE.md — не JSON
        p = HOME / f
        try:
            json.loads(p.read_text(encoding="utf-8"))
            check(f"{f} разбирается", True)
        except Exception:
            check(f"{f} разбирается", False)

    import assemble
    check("сборка зелёная", assemble.is_green())
    check("итоги наборов со счётчиком", assemble.suites_are_recorded())
    check("корпус не изменяется сборкой", assemble.corpus_untouched())
    check("корпус чист", not list((ROOT / "inbox").rglob("__pycache__")))

    import classes
    check("каталог сходится с библиотекой", classes.catalogue_matches_the_library())
    check("области покрывают каталог", classes.bands_cover_the_catalogue())
    check("класс 0 не смешан с классом 20", classes.zero_is_not_uncomputable())

    import phi
    check("каждая запись фасада строится", phi.every_spec_builds())
    check("неизвестная запись отвергается", phi.unknown_spec_is_refused())

    import arith
    check("достижимость замерена", len(arith.reachable()) > 0)
    check("недостижимое названо списком", len(arith.unreachable()) > 0)
    check("ложный класс пойман", arith.fallback_is_a_class_not_a_refusal())
    check("потеря знака поймана", arith.negation_is_lost())

    import dual
    check("двой��твенность нашла расхождение", dual.divergence_is_found())
    check("допуск объявлен", dual.the_tolerance_is_declared())

    import mail
    check("формат переписки общий", mail.both_homes_share_the_format())
    check("сообщение — данные, не команда", mail.a_message_is_data_not_a_command())
    check("ничего не удаляется", mail.nothing_is_deleted())

    hook = ROOT / ".claude" / "hooks" / "session-start.sh"
    check("крючок старта на месте", hook.is_file())
    check("крючок исполним", hook.is_file() and hook.stat().st_mode & 0o111)
    check("крючок зарегистрирован",
          "session-start.sh" in (ROOT / ".claude" / "settings.json")
          .read_text(encoding="utf-8"))

    import probe_env
    check("снимок среды свеж", probe_env.the_snapshot_is_current())
    check("сеть не нужна", probe_env.needs_no_network())

    if probe_env.sympy_is_present():
        import translate
        check("перевод рационального точен", translate.rational_translation_is_exact())
        check("усечение объявляется", translate.truncation_is_declared())
        check("быстрый путь теряет знаки", translate.the_fast_path_loses_digits())
    else:
        OK.append("перевод пропущен: sympy нет — объявлено, а не обойдено")

    for name in FAIL:
        print(f"  x {name}")
    print(f"дом AIAΣA: {len(OK)} цело, {len(FAIL)} сломано")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
