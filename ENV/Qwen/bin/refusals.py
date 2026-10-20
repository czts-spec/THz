#!/usr/bin/env python3
"""refusals — реестр названных отказов и крючок их явного пересмотра.

Слово автора 03.09: «названные отказы ПОДЛЕЖАТ ЯВНОМУ (рекурсивному)
ПЕРЕСМОТРУ после преодоления очередной стены… мы не можем себе позволить
закрыть 99 % математики потому, что мы её не знали». Устройство, не
дисциплина (routine-to-hook):

* `collect()` собирает все названные отказы библиотеки phicore: места
  `raise Unresolved / ExactRepresentationNotAvailable` в коде (по тексту
  подсказки, группами), отказные и пустые ячейки ФIOS (группами по
  семействам ключей), список «чего здесь пока нет» README. У каждой
  группы — устойчивый id (sha1 текста).
* `--check` (в воротах): сравнивает с последним снимком
  `state/refusals.json`. Исчезнувший отказ = ПРЕОДОЛЁННАЯ СТЕНА; тогда
  требуется ЯВНЫЙ ПЕРЕСМОТР всех оставшихся: генерируется заготовка
  `notes/пересмотр отказов после стены — <hash>. <дата>.md` с вердиктом
  «?» у каждой группы, и ворота красны, пока каждое «?» не заменено
  вердиктом. Пересмотр, снявший ещё отказы, вызовет следующий
  пересмотр — рекурсия устройством, до неподвижной точки.
* Новые названные отказы допускаются (они и есть честность), получают
  вердикт при следующем пересмотре; `--review` вызывает пересмотр
  принудительно; `--status` печатает счёт.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import pathlib
import re
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

NOT_MEASURED = (
    "реестр собирает отказы по тексту raise в phicore/ и по таблице ФIOS; "
    "отказы, выраженные иначе (return None, NotImplemented, TypeError), "
    "в него не входят — их доля не замерена",
    "вердикт пересмотра — суждение исполнителя, записанное явно; его "
    "верность устройство не проверяет, оно проверяет лишь его наличие",
)

sys.path.insert(0, str(HERE))
try:                                    # мерка «кого стена касается»
    import refusals_wall as WALL
except ImportError:                     # без мерки заготовка соберётся как прежде
    WALL = None

STATE = HERE.parent / "state" / "refusals.json"
NOTES = REPO / "notes"
PHICORE = REPO / "phicore"
VERDICT_RE = re.compile(r"\[([0-9a-f]{10})\].*?вердикт:\s*(.+?)\s*$")


# ---------------------------------------------------------------- сбор

_LIT = re.compile(r"(?:f|rf|fr)?(?:\"((?:[^\"\\\\]|\\\\.)*)\"|'((?:[^'\\\\]|\\\\.)*)')")


def _string_text(src: str) -> str:
    """Склейка строковых литералов внутри куска исходника."""
    parts = _LIT.findall(src)
    return re.sub(r"\s+", " ", "".join(a or b for a, b in parts)).strip()


def _str_consts(tree, src: str) -> dict:
    """Строковые константы модуля и классов: {имя: кусок исходника}.

    Имя класса входит в ключ (`ArcTrigValue._SUB_HINT`), чтобы одноимённые
    подсказки разных классов не путались. Берутся присваивания одному имени,
    в значении которых есть строковый литерал: этого хватает и на простую
    строку, и на склейку в скобках, и на f-строку."""
    out = {}

    def scan(body, prefix):
        for n in body:
            if (isinstance(n, ast.Assign) and len(n.targets) == 1
                    and isinstance(n.targets[0], ast.Name)):
                seg = ast.get_source_segment(src, n.value) or ""
                if _LIT.search(seg):
                    out[prefix + n.targets[0].id] = seg
            if isinstance(n, ast.ClassDef):
                scan(n.body, n.name + ".")

    scan(tree.body, "")
    return out


def _owners(tree) -> dict:
    """{узел raise: (имя класса, имя функции)} — место отказа в тексте."""
    out = {}
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            for m in ast.walk(n):
                if isinstance(m, ast.Raise):
                    cls, fn = out.get(id(m), ("", ""))
                    if isinstance(n, ast.ClassDef):
                        out[id(m)] = (n.name, fn)
                    elif not fn:
                        out[id(m)] = (cls, n.name)
    return out


def _const_key(a, cls: str = "") -> str | None:
    """Имя строковой константы, стоящей за аргументом, либо None.

    Разбирает три записи подсказки, и только их:

    * `ИМЯ` — константа модуля (`_POW_HINT`);
    * `self.ИМЯ` / `Объект.ИМЯ` — константа класса (`_SUB_HINT`);
    * `ИМЯ.format(…)` и `ИМЯ.join(…)` — ШАБЛОН, розданный по месту
      (`_ROOT_HINT.format(q=q)`, `hyper.py:729`).
    """
    if isinstance(a, ast.Name):
        return a.id
    if isinstance(a, ast.Attribute) and isinstance(a.value, ast.Name):
        return (cls + "." + a.attr) if a.value.id == "self" else a.attr
    if (isinstance(a, ast.Call) and isinstance(a.func, ast.Attribute)
            and a.func.attr in ("format", "join")):
        return _const_key(a.func.value, cls)
    return None


def _call_text(call, src: str, consts: dict, where: str, cls: str = "") -> str:
    """Текст отказа из вызова `raise …(…)`: ИМЯ РАЗРЕШАЕТСЯ ДО ЗНАЧЕНИЯ."""
    parts, unresolved = [], False
    for a in list(call.args) + [k.value for k in call.keywords]:
        seg = ast.get_source_segment(src, a) or ""
        if _LIT.search(seg):
            parts.append(_string_text(seg))
            continue
        key = _const_key(a, cls)
        if key and key in consts:
            parts.append(_string_text(consts[key]))
        elif isinstance(a, (ast.Name, ast.Attribute, ast.Call)):
            unresolved = True
    text = re.sub(r"\s+", " ", "".join(parts)).strip()
    if not text and unresolved:
        return (f"(подсказка именем, значение времени прогона; "
                f"отказ назван в {where})")
    return text


def _raise_sites():
    out = {}
    for p in sorted(PHICORE.glob("*.py")):
        s = p.read_text(encoding="utf-8")
        try:
            tree = ast.parse(s)
        except SyntaxError:
            continue
        consts, owners = _str_consts(tree, s), _owners(tree)
        for n in ast.walk(tree):
            if not (isinstance(n, ast.Raise) and isinstance(n.exc, ast.Call)):
                continue
            f = n.exc.func
            exc = f.id if isinstance(f, ast.Name) else getattr(f, "attr", None)
            if exc not in ("Unresolved", "ExactRepresentationNotAvailable"):
                continue
            cls, fn = owners.get(id(n), ("", ""))
            where = f"{p.name}:{cls or fn or '<модуль>'}"
            if cls and fn:
                where = f"{p.name}:{cls}.{fn}"
            text = _call_text(n.exc, s, consts, where, cls)
            out.setdefault((exc, text), []).append(f"{p.name}:{n.lineno}")
    items = []
    for (exc, text), sites in sorted(out.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        items.append({"kind": "код:" + ("отказ" if exc == "Unresolved" else "не-сведено"),
                      "text": text or "(без текста)",
                      "sites": sorted(sites, key=lambda s: (s.rsplit(":", 1)[0],
                                                            int(s.rsplit(":", 1)[1]))})
    return items


def _key_family(key) -> str:
    contract, cls = key
    if contract == "ALPHA":
        return "ALPHA"
    if contract == "STRUCT":
        return "STRUCT"
    if contract == "UNCOMP":
        return "класс 20" if cls == 20 else "класс 23"
    if cls == 0:
        return "класс 0"
    if cls in (21, 22):
        return "HYPER 21/22"
    if cls == 4:
        return "класс 4 (не веществен)"
    if cls == 18:
        return "класс 18 (сертификат по точке)"
    return "линейки 1–19"


def _fios_items():
    from phicore import Фios as F
    refuse, empty = {}, {}
    for a in F.LIVE:
        for b in F.LIVE:
            for op in "+-*/":
                cell = F.TABLE.get((a, b, op))
                if cell is None:
                    fam = (_key_family(a), _key_family(b))
                    empty.setdefault(fam, []).append(f"{a[0]}{a[1]}×{b[0]}{b[1]}«{op}»")
                elif cell[0] == "refuse":
                    refuse.setdefault(cell[1], []).append(f"{a[0]}{a[1]}×{b[0]}{b[1]}«{op}»")
    items = [{"kind": "ФIOS:отказ", "text": h, "sites": sorted(s)} for h, s in sorted(refuse.items())]
    items += [{"kind": "ФIOS:пусто", "text": f"{fa} × {fb}: ячейки пусты (иммунный отказ)",
               "sites": sorted(s)} for (fa, fb), s in sorted(empty.items())]
    return items


def _readme_items():
    s = (PHICORE / "README.md").read_text(encoding="utf-8")
    if "## Чего здесь пока нет" not in s:
        return []
    sec = s.split("## Чего здесь пока нет", 1)[1]
    items = []
    for m in re.finditer(r"^\* (.+?)(?=^\* |\Z)", sec, re.M | re.S):
        text = re.sub(r"\s+", " ", m.group(1)).strip()
        items.append({"kind": "README:нет", "text": text, "sites": ["phicore/README.md"]})
    return items


def collect() -> list:
    items = _raise_sites() + _fios_items() + _readme_items()
    for it in items:
        it["id"] = hashlib.sha1((it["kind"] + "|" + it["text"]).encode("utf-8")).hexdigest()[:10]
    return items


def set_hash(ids) -> str:
    return hashlib.sha1("\n".join(sorted(ids)).encode("utf-8")).hexdigest()[:8]


# ---------------------------------------------------------------- снимок


def load_state(path=STATE) -> dict:
    if not pathlib.Path(path).exists():
        return {"hash": None, "items": {}, "verdicts": {}, "reviewed": None,
                "base": None}
    return json.loads(pathlib.Path(path).read_text(encoding="utf-8"))


#: «оставить как было» для поля `base` — отличимо от None («основания нет»).
_KEEP = object()


def head_commit(repo=REPO) -> str | None:
    """Коммит, на котором дом стоит сейчас, — основание БУДУЩЕГО круга."""
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(repo),
                           capture_output=True, text=True)
    except OSError:
        return None
    sha = r.stdout.strip()
    return sha if r.returncode == 0 and len(sha) == 40 else None


def save_state(items, verdicts, reviewed, path=STATE, base=_KEEP) -> None:
    ids = [it["id"] for it in items]
    if base is _KEEP:
        base = load_state(path).get("base")
    data = {"hash": set_hash(ids),
            "base": base,
            "items": {it["id"]: {"kind": it["kind"], "text": it["text"], "sites": it["sites"]}
                      for it in items},
            "verdicts": {i: v for i, v in verdicts.items() if i in ids},
            "reviewed": reviewed}
    pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=1) + "\n",
                                  encoding="utf-8")


# ---------------------------------------------------------------- пересмотр


def review_path(h: str, notes=NOTES) -> pathlib.Path | None:
    for p in sorted(pathlib.Path(notes).glob(f"пересмотр отказов после стены — {h}.*.md")):
        return p
    return None


def _splitter(state=STATE):
    if WALL is None:
        return None
    try:
        s = WALL.Split(base=WALL.snapshot_base(REPO, state))
    except Exception:
        return None
    return s if s.usable else None


def skeleton(items, removed: dict, prev_verdicts: dict, h: str, notes=NOTES,
             verbatim: dict | None = None, split=None) -> pathlib.Path:
    day = time.strftime("%Y-%m-%d")
    p = pathlib.Path(notes) / f"пересмотр отказов после стены — {h}. {day}.md"
    if split is None:
        split = _splitter()
    wall = WALL.wall_line(removed) if WALL else ""
    stamp = f"{h}. {day}"
    verbatim = verbatim or {}
    body, carried, hand, kept = [], 0, 0, 0
    for it in items:
        prev = prev_verdicts.get(it["id"])
        hint = f" (прежний вердикт: {quoted(prev)})" if prev else ""
        body.append(f"- [{it['id']}] {it['kind']} · „{it['text']}“{hint}")
        body.append(f"  места ({len(it['sites'])}): " + ", ".join(it["sites"][:12])
                    + (" …" if len(it["sites"]) > 12 else ""))
        if it["id"] in verbatim:
            body.append(f"  вердикт: {verbatim[it['id']]}")
            kept += 1
            continue
        v = (WALL.carried_verdict(it["sites"], prev, wall, stamp, split)
             if split is not None and prev else None)
        if v:
            body.append(f"  вердикт: {v}")
            carried += 1
        else:
            body.append("  вердикт: ?")
            hand += 1
    lines = [f"# Пересмотр отказов после стены — {h}", "",
             "Слово автора 03.09: названные отказы подлежат явному (рекурсивному)",
             "пересмотру после преодоления очередной стены. Заготовка собрана",
             "устройством `ENV/AIAΣA/bin/refusals.py`; ворота красны, пока у какой-либо",
             "группы вердикт не проставлен. Вердикт — одна строка: «стоит: <почему>»,",
             "«снят: <чем>» или «сужен: <до чего>». Провенанс — [вычислено здесь].", ""]
    if split is not None:
        lines += [
            f"Разделение сделано устройством `ENV/AIAΣA/bin/refusals_wall.py` "
            f"замером правок круга (`git status`, `git diff HEAD`):",
            f"перенесено прежних вердиктов — **{carried}**, сохранено дословно "
            f"из прошлой сборки этого круга — **{kept}**, оставлено на суждение "
            f"— **{hand}**.",
            "Перенос идёт ТОЛЬКО туда, где место отказа этим кругом не правилось;",
            "где правлена сама строка отказа либо несущая его функция — стоит «?»,",
            "и ворота красны, пока там не сказано «снят/сужен/стоит» с адресом.", ""]
    else:
        lines += ["Правок круга не видно (`git diff HEAD` пуст либо мерка "
                   "недоступна): перенос НЕ делался, все группы на суждение — "
                   "это факт о мерке, не о доме.", ""]
    lines += ["## Преодолённая стена — отказы, исчезнувшие из реестра", ""]
    if removed:
        for i, it in removed.items():
            lines.append(f"- [{i}] {it['kind']} · „{it['text']}“")
    else:
        lines.append("- (первый полный пересмотр: снимка прежде не было)")
    lines += ["", f"## Оставшиеся названные отказы — {len(items)} групп", ""]
    lines += body
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


_HEADS = ("снят:", "сужен:", "стоит:")


def quoted(v: str) -> str:
    v = v.strip()
    for h in _HEADS:
        if v.startswith(h):
            rest = v[len(h):].strip()
            return v if rest.startswith("„") else f"{h} „{rest}“"
    return v if v.startswith("„") else f"„{v}“"


def parse_verdicts(p: pathlib.Path) -> dict:
    out = {}
    for ln in p.read_text(encoding="utf-8").splitlines():
        m = re.match(r"^- \[([0-9a-f]{10})\]", ln)
        if m:
            cur = m.group(1)
            continue
        m = re.match(r"^\s+вердикт:\s*(.*?)\s*$", ln)
        if m and "cur" in dir():
            out[cur] = m.group(1)
    return out

_REVIEW_NAME = re.compile(r"пересмотр отказов после стены — ([0-9a-f]{6,10})\.\s*"
                         r"(\d{4}-\d{2}-\d{2})\.md$")


def review_order(p: pathlib.Path) -> tuple:
    m = _REVIEW_NAME.search(p.name)
    return (m.group(2) if m else "", p.name)


def harvest_fresh(notes=NOTES, state=STATE) -> dict:
    known = load_state(state).get("verdicts", {}) or {}
    out, absorbed = {}, set()
    for p in sorted(pathlib.Path(notes).glob("пересмотр отказов после стены — *.md"), key=review_order):
        for i, v in parse_verdicts(p).items():
            if v.strip() in ("", "?"):
                continue
            if known.get(i) == v:
                absorbed.add(i)
                out.pop(i, None)
                continue
            if i in absorbed:
                continue
            out[i] = v
    return out


def check(force: bool = False, state=STATE, notes=NOTES) -> int:
    items = collect()
    ids = {it["id"] for it in items}
    prev = load_state(state)
    prev_ids = set(prev.get("items", {}))
    removed = {i: prev["items"][i] for i in prev_ids - ids}
    h = set_hash(ids)
    verdicts = dict(prev.get("verdicts", {}))
    due = force or bool(removed) or prev.get("hash") is None
    if not due:
        if prev.get("hash") != h:
            save_state(items, verdicts, prev.get("reviewed"), state)
        pending = [it for it in items if it["id"] not in verdicts]
        if pending:
            p_now = review_path(h, notes)
            if p_now is not None:
                got = parse_verdicts(p_now)
                fresh = {i["id"]: got[i["id"]] for i in pending
                         if got.get(i["id"], "?").strip() not in ("", "?")}
                if fresh:
                    verdicts.update(fresh)
                    save_state(items, verdicts, p_now.name, state)
                    pending = [it for it in items if it["id"] not in verdicts]
                    print(f"принято вердиктов из пересмотра {p_now.name}: {len(fresh)}")
        print(f"отказов: {len(items)} групп; без вердикта (новые): {len(pending)}; стена не преодолевалась")
        return 0
    p = review_path(h, notes)
    if p is None:
        fresh = harvest_fresh(notes, state)
        p = skeleton(items, removed, {**verdicts, **fresh}, h, notes,
                     verbatim=fresh, split=_splitter(state))
        print(f"СТЕНА ПРЕОДОЛЕНА: исчезло отказов {len(removed)}; требуется явный пересмотр "
              f"{len(items)} групп — заготовка: {p.name}")
        return 1
    split = _splitter(state)
    if split is not None and WALL is not None:
        res = WALL.run(p, WALL.wall_line(removed), {}, split,
                       sites_by_id={it["id"]: it["sites"] for it in items})
        if res["count"]["вне"] or res["count"]["файл"]:
            print(f"перенесено устройством: вне правки круга "
                  f"{res['count']['вне']}, место не тронуто "
                  f"{res['count']['файл']}; на суждение {len(res['hand'])}")
    fresh = harvest_fresh(notes, state)
    if fresh:
        try:
            import refusals_merge as MERGE
            got = MERGE.merge(p, fresh)
            if got["filled"]:
                print(f"перенесено из соседних файлов пересмотра: "
                      f"{len(got['filled'])}")
        except ImportError:
            pass
    got = parse_verdicts(p)
    missing = [it["id"] for it in items if got.get(it["id"], "?").strip() in ("", "?")]
    if missing:
        print(f"пересмотр не завершён: без вердикта {len(missing)} из {len(items)} — {p.name}")
        return 1
    verdicts.update({i: got[i] for i in ids})
    save_state(items, verdicts, p.name, state, base=head_commit())
    print(f"пересмотр завершён: {len(items)} вердиктов, снимок {h} записан ({p.name})")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--review", action="store_true", help="принудительный пересмотр")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for it in collect():
            print(f"[{it['id']}] {it['kind']} ({len(it['sites'])}) · {it['text'][:110]}")
        return 0
    if a.status:
        items = collect(); prev = load_state()
        v = prev.get("verdicts", {})
        print(f"названных отказов: {len(items)} групп; с вердиктом: "
              f"{sum(1 for it in items if it['id'] in v)}; снимок: {prev.get('hash')}; "
              f"последний пересмотр: {prev.get('reviewed')}")
        return 0
    return check(force=a.review)


if __name__ == "__main__":
    raise SystemExit(main())
