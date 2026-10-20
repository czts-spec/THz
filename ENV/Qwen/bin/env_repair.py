#!/usr/bin/env python3
"""Самовосстановление среды при старте свежего контейнера.

Свежий контейнер приходит с двумя изъянами, из-за которых наборы среды
краснеют, хотя репозиторий цел (найдено 02.09 на воротах
total_nonFraction):

1. усечённый клон гита — наборы истории (`bb_archive`,
   `transcript_gap`) не видят коммитов;
2. пустые срезы стенограмм во временном каталоге — `line_names` и
   родня меряют по ним живым проходом.

Инструмент чинит оба изъяна и молчит, когда чинить нечего. Правило
дома держится: код возврата 0 достижим без сети — сеть нужна только
самому акту дочитки истории, и если её нет, остаток называется
поимённо, а не замалчивается.

Дом: своё хранится в своём доме (слово автора 02.09); крючок старта
`.claude/hooks/session-start.sh` лишь вызывает этот инструмент.
"""
import pathlib
import subprocess
import sys

NOT_MEASURED = (
    "outbox-каталоги домов не проверяются: пустой каталог держится "
    "в гите файлом .gitkeep, это починено в самом репозитории",
    "здоровье наборов не прогоняется отсюда: это работа ворот "
    "proto/gate.sh; инструмент лишь возвращает среде исходное",
)

ROOT = pathlib.Path(__file__).resolve().parents[3]
HERE = pathlib.Path(__file__).resolve().parent


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(ROOT), *args],
                          capture_output=True, text=True)


def clone_is_shallow() -> bool:
    r = _git("rev-parse", "--is-shallow-repository")
    return r.stdout.strip() == "true"


def unshallow() -> bool:
    """Дочитать историю; правда возврата — состояние, не намерение."""
    _git("fetch", "--unshallow", "origin")
    return not clone_is_shallow()


def slices_are_whole() -> bool:
    sys.path.insert(0, str(HERE))
    import deploy_slices
    n, b = deploy_slices.acceptance()
    return n == 181 and b == 40_887_618


def deploy() -> bool:
    r = subprocess.run([sys.executable, "-B",
                        str(HERE / "deploy_slices.py")],
                       capture_output=True, text=True)
    return r.returncode == 0


if __name__ == "__main__":
    remainder = []

    if clone_is_shallow():
        if unshallow():
            print("клон дочитан до полной истории")
        else:
            remainder.append("клон остался усечённым (нет сети?) — "
                             "наборы истории будут красны")
    else:
        print("история полна")

    if slices_are_whole():
        print("срезы стенограмм на месте (181 файл)")
    elif deploy():
        print("срезы стенограмм развёрнуты заново (181 файл, "
              "40 887 618 байт)")
    else:
        remainder.append("срезы не развернулись — "
                         "python3 -B ENV/AIAΣA/bin/deploy_slices.py")

    # Сбой площадки: зеркала среды нет — вернуть из резервной копии репозитория
    # (слово автора 25.09, c7c19b5c42), реестр объединяется, не затирается.
    _mirror = pathlib.Path.home() / ".claude" / "aiasa" / "state" / "course.json"
    if not _mirror.exists():
        _r = subprocess.run([sys.executable, "-B", str(pathlib.Path(__file__).resolve().parent / "env_backup.py"), "--restore"],
                            capture_output=True, text=True)
        print("зеркало среды из резервной копии: " + (_r.stdout.strip() or _r.stderr.strip()[-120:]))

    # Правила дома на уровне среды (слово автора 25.09, e99628ba0f): ~/.claude —
    # CLAUDE.md, крючки, зеркало реестра; идемпотентно, чужое не трогает.
    _inst = subprocess.run([sys.executable, "-B", str(pathlib.Path(__file__).resolve().parent / "install_env.py"), "--quiet"],
                            capture_output=True, text=True)
    print("правила дома на уровне среды: ~/.claude — " + ("установлены" if _inst.returncode == 0 else "НЕ установлены: " + _inst.stderr.strip()[-120:]))

    for r in remainder:
        print("ОСТАТОК:", r)
    sys.exit(1 if remainder else 0)
