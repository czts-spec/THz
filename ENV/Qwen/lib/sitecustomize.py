# sitecustomize — среда Qwen (адаптация ФIOS-обвязки из InvesePolar).
# Контейнерная среда не даёт OpenBLAS прочитать размер кэша L2 -> ложный
# warning "could not determine the L2 cache size" при ПЕРВОМ импорте numpy.
# Вычислений это не касается: OpenBLAS просто берёт значение по умолчанию.
# Заглушка точечная: глушим fd 2 только на время импорта numpy, восстанавливаем
# немедленно — настоящие ошибки и трейсбеки видны полностью (канон честности).
import os, sys

def _mute_first_import():
    try:
        import numpy  # noqa: F401
    except Exception:
        pass

if "numpy" not in sys.modules:
    _r, _w = os.pipe()
    _old = os.dup(2)
    os.dup2(_w, 2)
    try:
        _mute_first_import()
    finally:
        os.dup2(_old, 2)
        os.close(_w); os.close(_r); os.close(_old)
