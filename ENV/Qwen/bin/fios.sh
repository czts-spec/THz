#!/bin/sh
# ФIOS-обвязка для среды Qwen (носитель: InvesePolar/phicore, скопирован только для чтения).
# 1) подавляет ложный OpenBLAS-шум контейнера (L2 cache недоступен в sandbox);
#    OpenBLAS пишет warning напрямую в fd 2 — перекрываем devnull на время импорта;
# 2) добавляет точную арифметику phicore (Q, AlgebraicNumber, ...) в PYTHONPATH.
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
_FIOS_LIB="$(CDPATH= cd -- "$(dirname -- "$0")/../lib" 2>/dev/null && pwd)"
[ -n "$_FIOS_LIB" ] && export PYTHONPATH="$_FIOS_LIB${PYTHONPATH:+:$PYTHONPATH}"
# Пред-грев: импорт numpy в отдельном процессе, где stderr -> /dev/null.
# OpenBLAS пишет warning на уровне fd 2 (не через sys.stderr), поэтому глушим
# именно файловый дескриптор. Реальный запуск идёт с НЕТРОНУТЫМ stderr:
# настоящие ошибки и трейсбеки видны полностью (канон: не проглатываем суть).
python3 -c "import numpy" >/dev/null 2>/dev/null || true
if [ -n "$1" ]; then
  exec python3 "$@"
fi
