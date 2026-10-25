# ФIOS-обвязка для среды Qwen (адаптация из InvesePolar, ТОЛЬКО ДЛЯ ЧТЕНИЯ источника)

Назначение: OpenBLAS в контейнерной среде выдаёт ложный шум
`WARNING - could not determine the L2 cache size on this system, assuming 256k`.
Проверено вживую: warning появляется при первом импорте numpy на уровне fd 2,
результаты вычислений не искажаются (OpenBLAS берёт дефолт 256k). Шум маскирует
настоящие сообщения об ошибках — по канону убран.

Состав:
- `lib/phicore/` — копия носителя истины InvesePolar/phicore (Фios.py и т.д.),
  только для чтения; даёт точную арифметику (Q, AlgebraicNumber, ...) через PYTHONPATH.
  Двойное имя модуля сохранено как у автора: phicore.fios is phicore.Фios.
- `lib/sitecustomize.py` — точечная заглушка: глушит fd 2 только на время первого
  импорта numpy, восстанавливает немедленно. Настоящие трейсбеки видны полностью.
- `bin/fios.sh` — обвязка запуска: экспорт тред-переменных, PYTHONPATH=lib, exec python3.

Применение в среде Qwen:
    sh ENV/Qwen/bin/fios.sh скрипт.py [аргументы]
или вручную:
    export PYTHONPATH=$PWD/ENV/Qwen/lib

Источник скопирован из https://github.com/czts-spec/InvesePolar (phicore/) —
репозиторий только для чтения, изменения туда не вносятся.
Дата адаптации: 2026-10-25. Прогон проверки: numpy sum ok, ФIOS Q(1,3)+Q(1,6)==Q(1,2) True,
ZeroDivisionError показан целиком (заглушка его не съела).
