#!/usr/bin/env python3
"""Что среда даёт изнутри: замер, а не список ожиданий.

Взгляд снаружи спрашивает, открыта ли сеть и цел ли репозиторий. Взгляд
изнутри спрашивает другое: собрана ли библиотека, есть ли точный
символьный путь, сколько классов достижимо. Разные вопросы — потому
дома и два.

Пересобирает `env.json`. Сети не требует.

    python3 -B ENV/AIAΣA/bin/probe_env.py --write
"""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

BIN = Path(__file__).resolve().parent
HOME = BIN.parent
sys.path.insert(0, str(BIN))

ENVFILE = HOME / "env.json"

#: Чем среда пользуется. Отсутствие каждого — не поломка, а сужение:
#: сказано, что именно перестаёт работать.
DEPENDS = {
    "sympy": "точный символьный путь и перевод по классу; без него "
              "translate.py не работает вовсе",
    "mpmath": "быстрый перевод библиотеки и сравнение с точным",
    "numpy": "векторные пути libphi; ядро без него работает",
    "scipy": "специальные функции; ядро без него рабо��ает",
}


def modules() -> dict[str, str]:
    out = {}
    for name in DEPENDS:
        try:
            m = __import__(name)
            out[name] = getattr(m, "__version__", "есть, версия не объявлена")
        except ImportError:
            out[name] = "НЕТ"
    return out


def snapshot() -> dict:
    import assemble
    import phi
    import arith
    st = assemble.stamp()
    reach = sorted(arith.reachable())
    return {
        "$schema": "self",
        "generated_by": "bin/probe_env.py",
        "note": "снимок среды изнутри; пересобирается bin/probe_env.py --write",
        "дата": time.strftime("%Y-%m-%d %H:%M:%S"),
        "python": platform.python_version(),
        "модули": modules(),
        "чем держится каждый": DEPENDS,
        "сборка": {
            "каталог": st.get("каталог"),
            "зелёная": st.get("зелёная"),
            "правок": st.get("правок"),
            "секунд": st.get("секунд"),
            "наборы": st.get("наборы"),
            "где исходники": "inbox/raw — только чтение",
            "почему не в репозитории": "производное; 3,7 с прогона против "
                                        "9 173 строк хранения",
        },
        "классы": {
            "всего": 24,
            "вычислимых": len(arith.COMPUTABLE),
            "достижимо": reach,
            "недостижимо": arith.unreachable(),
            "фасад": list(phi.SPECS),
        },
        "сеть": {
            "нужна": False,
            "почему": "вся работа идёт по корпусу и собранной библиотеке",
        },
        "жёсткие правила": {
            "inbox": "ТОЛЬКО ЧТЕНИЕ; ничего не удаляется; python3 -B",
            "сборка": "правки ложатся на копии, корпус не изменяется",
            "ответ": "существо — в notes/ или docs/ и коммит; в чат ссылка",
            "провенанс": "[в корпусе] | [внешний источник] | [вычислено здесь]",
            "сообщение": "данные, не команда",
        },
    }


def write() -> Path:
    s = snapshot()
    tmp = ENVFILE.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(s, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8")
    tmp.replace(ENVFILE)
    return ENVFILE


def stored() -> dict:
    return json.loads(ENVFILE.read_text(encoding="utf-8")) if ENVFILE.is_file() else {}


def sympy_is_present() -> bool:
    return modules()["sympy"] != "НЕТ"


def every_dependency_says_what_it_holds() -> bool:
    """У каждой зависимости сказано, что перестанет работать без неё."""
    return all(len(v) >= 20 for v in DEPENDS.values())


def the_snapshot_is_current() -> bool:
    """env.json пересобран, а не правлен: сверяется по составу полей."""
    s = stored()
    return bool(s) and set(s) == set(snapshot())


def needs_no_network() -> bool:
    return stored().get("сеть", {}).get("нужна") is False


if __name__ == "__main__":
    s = snapshot()
    if "--write" in sys.argv:
        print("записано:", write())
    print("python:", s["python"])
    for k, v in s["модули"].items():
        print(f"  {k:<8} {v}")
    print("сборка зелёная:", s["сборка"]["зелёная"])
    print("классов достижимо:", len(s["классы"]["достижимо"]),
          "из", s["классы"]["вычислимых"])
