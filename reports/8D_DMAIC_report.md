# 8D / DMAIC Report: Etch Module Excursions and Pump Failure

**Area:** Plasma etch module, chambers ETCH01-ETCH04 (simulated)
**Period:** days 0-60
**Author:** Rachana Santhebachahalli Umesha

| 8D step | DMAIC phase | Summary |
| --- | --- | --- |
| D1 Team | Define | Process, equipment, metrology and maintenance engineering |
| D2 Problem | Define | Etch rate went out of control on two chambers; one chamber lost 36 h to an unplanned pump trip; about 3% of wafers fail non-uniformity at the old recipe |
| D3 Containment | Measure | Hold lots from the alarmed chambers for extra metrology; move WIP to ETCH01 |
| D4 Root cause | Analyse | Pressure offset after PM (ETCH02); chiller degradation (ETCH03); accelerating pump bearing wear (ETCH04); recipe sitting where the RF × temperature interaction is unfavourable |
| D5 Corrective actions | Improve | Valve re-calibration; chiller repair; predictive pump swap; DOE recipe |
| D6 Verification | Improve | Confirmation runs; SPC back in control; forecast checked against the actual failure |
| D7 Prevention | Control | PM checklist and post-PM qualification; daily wear-model refit; SPC on sensors |
| D8 Closure | Control | FMEA RPN 1035 → 351; lessons shared |

## D2 / Define: problem statement

- **ETCH02:** from day 32 the etch rate ran 3.2 nm/min high, and Cpk fell from 1.78 to 1.31, below the 1.33 target.
- **ETCH03:** from day 20 the etch rate drifted low. Cpk fell from 1.71 to 1.03, and only 59.8% of wafers met both specs in the last 10 days.
- **ETCH04:** the dry pump tripped on high vibration on day 51.0, taking the chamber down for 36 hours.
- **All chambers:** at the process of record, about 3% of wafers exceed the 4% non-uniformity limit (Ppk 0.64).

## D4 / Measure and Analyse: root cause

1. **Detection.** The X-bar chart flagged ETCH02 within 0.5 days and ETCH03 on day 24.5, with continuous alarms from day 32. ANOVA confirmed the chambers no longer matched (days 50-60: F = 453, p < 0.001). Tukey HSD showed ETCH02 and ETCH03 differ from the others, while ETCH01 and ETCH04 do not (p = 0.98).
2. **Sensor trace.** Daily sensor means isolated the cause on each chamber. ETCH02 pressure stepped up 1.6 mTorr at the day-32 PM. ETCH03 chuck temperature climbed about 0.25 °C per day from day 20. Other sensors were unchanged.
3. **5 Whys, ETCH02.** Etch rate high → pressure high → throttle valve offset → valve not re-zeroed at PM → no calibration step or post-PM check in the PM procedure (**root cause**).
4. **Pump.** Vibration rose with an accelerating trend. The ODE wear model fitted it with an RMSE of 0.12 mm/s, consistent with progressive bearing wear.
5. **Recipe.** The DOE showed non-uniformity is driven mainly by the RF power × chuck temperature interaction (t = 11.9). The old recipe sat where this interaction is unfavourable.

## D5-D6 / Improve: actions and verification

| Action | Verification |
| --- | --- |
| Re-calibrate the ETCH02 throttle valve | Pressure back to 30.0 mTorr; X-bar back inside limits |
| Repair the ETCH03 chiller | Chuck temperature back to 60 °C |
| Replace the ETCH04 pump on a schedule set by the wear model | The forecast made on day 40 predicted day 51.3 (90% range 51.0-51.7); actual day 51.0. A planned swap on day 49 avoids the unplanned trip |
| New recipe: 1428 W / 32 mTorr / 55 °C | 200 confirmation wafers: non-uniformity 3.39% → 2.82%; Ppk 0.64 → 1.30; 100% in spec |

## D7 / Control: prevent recurrence

- Add valve calibration and two qualification wafers to the post-PM checklist.
- Use three-way X-bar / R charts with Western Electric rules 1-4 on etch rate. Use control charts on chuck temperature and pressure, so a drift is caught at the sensor before it reaches the product.
- Refit the pump wear model daily, and plan a pump swap at the 5th-percentile forecast failure day minus 1 day.
- Run automated SQL extraction and a data-quality report on every run, replacing manual spreadsheet tracking.

## FMEA (top items; full table in outputs/fmea_ranked.csv)

| Failure mode | S | O | D | RPN | RPN after |
| --- | --- | --- | --- | --- | --- |
| Pump trips during processing | 9 | 5 | 6 | 270 | 36 |
| Sensor dropouts / duplicate records | 5 | 5 | 7 | 175 | 50 |
| Chuck temperature drift | 7 | 4 | 6 | 168 | 84 |
| Pressure offset after PM | 7 | 4 | 5 | 140 | 28 |
| Non-uniformity above spec | 6 | 6 | 3 | 108 | 36 |
