#!/usr/bin/env bash
# Обновление резервного зеркала workspace (.git -> ENV/Qwen/mirror/workspace.git)
# Происхождение скрипта: [вычислено здесь], сессия 2026-10-19.
set -euo pipefail

SRC="${1:-/workspace/.git}"
HOME_DIR="$(cd "$(dirname "$0")" && pwd)"
DST="$HOME_DIR/mirror/workspace.git"

if [ ! -d "$SRC" ]; then
    echo "источник не найден: $SRC" >&2
    exit 1
fi

rm -rf "$DST.tmp"
cp -a "$SRC" "$DST.tmp"
mv "$DST" "$DST.prev" 2>/dev/null || true
mv "$DST.tmp" "$DST"
[ -d "$DST.prev" ] && rm -rf "$DST.prev"

echo "зеркало обновлено: $DST"
git -C "$DST" fsck --no-progress >/dev/null && echo "целостность: OK"
