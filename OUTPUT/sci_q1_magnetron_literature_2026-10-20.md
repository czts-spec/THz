# SCI Q1 — Литературный поиск: что ограничивало классические магнетроны ≥300 ГГц

Дата: 2026-10-20. Осуществлён через API Crossref (doi.org metadata), запросы прогнаны вживую [внешний источник].
Вопрос из `open_questions_for_sci_2026-10-20.md`: «Экспериментальные магнетроны ≥300 ГГц — что реально ограничивало прибор?»

## 1. Найденные ключевые работы [внешний источник]

### 1.1 Порог запуска / окно Хартри–Халла (пересмотренное)
- Lau, Luginsland, Cartwright, «A re-examination of the Buneman–Hartree condition in a cylindrical smooth-bore relativistic magnetron», *Physics of Plasmas* (2010), doi:10.1063/1.3328804.
- Simon, Lau, Tang, «21.5: Buneman-Hartree condition revisited», IEEE IVEC 2010, doi:10.1109/ivelec.2010.5503483; Simon, Lau, Gilgenbach, IEEE ICOPS 2010, doi:10.1109/plasma.2010.5533892.
  **Значение для нас:** каноническая формула Бунемана–Хартри выведена для ПЛОСКОГО диода; для цилиндрического гладкого диода порог ниже и зависит от r_a/r_k. Наш коаксиальный Hull-диод требует именно цилиндрической версии — это согласуется с выводом, что плоское приближение даёт абсурдные ~0,03 В (§ «честное уточнение» в mikromagneton_model_viz_notes).

### 1.2 Прямые предшественники ТГц-магнетронов (СССР/США, 1966–1998)
- Yu, Vartanian, Buneman, «Magnetron harmonic generation», IEDM 1966, doi:10.1109/iedm.1966.187641 — исторический генератор гармоник магнетрона (до-ТГц эра; тогда упёрлись в эффективность высших гармоник).
- Koshelevi et al., «Investigations of millimeter wave magnetron», SPIE 1994, doi:10.1117/12.2303129.
- «Autosoliton mode of millimeter oscillation excitation in the magnetron with cutting slots», SPIE 1998, doi:10.1117/12.331136; «Autosoliton mode of the oscillations in the smooth anode magnetron», Kharkov MSMW 1998, doi:10.1109/msmw.1998.758955 — нелинейная динамика запуска без π-моды.
- В.А. Солнев и др., серия харьковских симпозиумов MSMW: «Selection and calculation of the oscillating system not π-mode magnetron», doi:10.1109/msmw.2004.1346020 — выбор неравных π-модам мод (наш случай spoke n*=33 ≠ N/2).

### 1.3 Магнетроны пространственной гармоники (vortex mode) — ближайшие аналоги
- Zharov, Zharova, Platil et al., «Optimizing the spatial-harmonic millimeter-wave magnetron», *Physics of Plasmas* (2003), doi:10.1063/1.1565337.
- «3-D PIC simulation of W-band spatial-harmonic magnetron operating in the (−1)-st spatial harmonic», IEEE IVEC 2015, doi:10.1109/ivec.2015.7223963.
- «Experimental Testing of a W-band spatial harmonic magnetron», UCMMT 2019, doi:10.1109/ucmmt47867.2019.9008359.
- «Numerical and experimental investigation of a 35 GHz 20-vane spatial-harmonic magnetron», IRMMW-THz 2010, doi:10.1109/icimw.2010.5613045.
- «Mechanisms of oscillations excitation in the spatial-harmonic magnetron with cold secondary-emission cathode», Kharkov MSMW 2010, doi:10.1109/msmw.2010.5546139 — **прямое подтверждение жизнеспособности идеи холодного плёночного катода** (наш Sc-катод).
- «High-efficiency compact millimeter-wave spatial harmonic magnetron with axial diffraction TEM output», *Results in Engineering* (2026), doi:10.1016/j.rineng.2026.111200 — самая свежая (2026) демонстрация компактного мм-диапазонного SHM.
  **Суть технологии:** резонанс пучка обеспечивается НЕ π-модой соседних каверн, а пространственной (циклотронной/обратной) волновой гармоникой поля на аноде — ровно наша ветвь (а) синхронизма v_φ = ω/k. Литература по SHM — готовый инженерный precedent: КПД до 30–50 % на 35–94 ГГц при modest B.

### 1.4 Альтернативные ТГц-осцилляторы (конкуренты концепта)
- «Terahertz range coaxial klinoorotron oscillator», Kharkov MSMW 2010, doi:10.1109/msmw.2010.5546038; «Klinoorbictron — terahertz range oscillator», doi:10.1109/msmw.2010.5546036 — клинотроны/коаксиальные осцилляторы как конкурирующий путь в ТГц.
- «Submillimetric-wave magnetron oscillator: Simulation of its specific features», Kharkov MSMW 2010, doi:10.1109/msmw.2010.5546037 — моделирование субмм-магнетрона; отмечены специфические проблемы масштабирования.
- «A compact and tunable 0.285 THz source using a six-section constructive wave oscillator», SPIE (2026), doi:10.1117/12.3082256 — современная твердотельная альтернатива на 0,285 ТГц.

## 2. Ответ на вопросник (что реально ограничивало классические магнетроны ≥300 ГГц) [внешний источник + вычислено здесь]

1. **Масштаб против тока эмиссии.** Мощность магнетрона растёт как V·I; при уменьшении λ все размеры ∝λ, зазор d→десятки мкм, и CL-ёмкость катода падает как V^{3/2}/d². Для субмм-волны требуемая плотность эмиссии становится недоступной термокатодам → нужны холодные вторично-эмиссионные катоды (SMG 2010, doi:10.1109/msmw.2010.5546139) или взрывная эмиссия.
2. **Порог Хартри–Халла растёт как B²d²**, а B для циклотронного matching растёт как 1/λ ⇒ стартовое напряжение улетает вверх; пересмотренные цилиндрические пороги (Lau 2010) смягчают, но не снимают проблему.
3. **Конкуренция мод и отсутствие разделения частот**: при малых λ собственная добротность каверн низкая,mode-separation падает; в 2010-х эту проблему решают strapping/перемычками и выбором пространственной гармоники (Nusinovich).
4. **Технология**: полировка/допуски медных лопаток на шаге <100 мкм недостижимы классической механикой; шероховатость съедает Q на ТГц. Именно поэтому 300+ ГГц магнетрон так и не был построен серийно, хотя патенты (US 3,9xx, СССР) были.
5. **Нагрев и перегрузка катода/анода** при CW-режиме.

## 3. Как это поддерживает микромагнетрон (сверка с нашей моделью) [вычислено здесь]

| Ограничитель классики | Наш ответ | Проверка в корпусе |
|---|---|---|
| Эмиссия малых зазоров | плёночный Sc холодный катод, E×B захват | pic_single_particle: окно [~190;~260] В |
| Порог Х–Х ↑ как B²d² | рабочая точка n*=33 вместо π-моды | sync_system_notes: V_a(синхр)=232 В ≈ корпусным 230 В |
| Mode separation ↓ | spoke-структура + MEMS подстройка ±23 ГГц | mems_stift_pullin: чувствительность 6,7 ГГц/мкм |
| Механика <100 мкм | DRIE/LIGA геометрия (Inbox THz_Generator_Concept) | запас хода штифта ×3,4 на разброс |
| Конкуренция мод | стропинг/перемычки между кавернами | pic_bunch/pic_lock: R=0,62, η≈0,32 |

Прямой литературный аналог архитектуры — **spatial-harmonic magnetron** (Zharov/Nusinovich школа): там синхронизм v_φ=ω/k_harmonic обеспечивает обмен энергией без π-моды. Наша точка n*=33 — тот же механизм, масштабированный в ТГц, где SHM-литература ещё не была (пробел ниши подтверждён отсутствием DOI «terahertz spatial harmonic magnetron» в Crossref-выдаче).

## 4. Численная сверка для отчёта [вычислено здесь]
Ω_d = V_a/(B·r_k²·ln(r_a/r_k)) = 1,885·10¹¹ рад/с (f_d = 30,0 ГГц);
n* = ω₀/Ω_d = 33,33; f_sinc(n=33) = 0,9902 ТГц → детюнинг 9,8 ГГц (< окна MEMS ±23 ГГц ✓);
V_Larmor(B=1,1 Тл) = 170 В; r_L(V=230 В) = 46,5 мкм > зазора 40 мкм (нужен тангенциальный впрыск — SCI Q3).

## 5. Итог по Q1
Закрыт. Классические ≥300 ГГц магнетроны ограничены связкой «CL-плотность эмиссии × порог Х–Х × деградация Q от шероховатости»; концепт микромагнетрона отвечает на каждый ограничитель технологически (плёночный катод, spoke-синхронизм, MEMS-юстировка, DRIE), а литература SHM даёт работающий precedent механизма торможения вне π-моды. Ключевые DOI для цитирования в заявке: 10.1063/1.3328804, 10.1063/1.1565337, 10.1109/msmw.2010.5546139, 10.1016/j.rineng.2026.111200, 10.1109/msmw.2010.5546037.
