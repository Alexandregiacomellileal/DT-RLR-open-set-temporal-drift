# DT-RLR repository update for release v1.1.0

This update adds the experiments introduced after peer-review feedback and used in the Measurement-targeted Revised 3 manuscript.

## Added analyses

- Safe-gate sensitivity: z_gate = 3, 4, 5.
- Slow-reference update sensitivity: eta = 0.10, 0.20, 0.40.
- Module ablation: drift trigger with slow adaptation but without the fast-context gate.
- Multiple simultaneous UNKNOWN fault types (2 and 3 held-out types).
- Runtime and memory-overhead audit for DT-RLR reference management.

## Main new result files

- `results/paper2_revised3_sensitivity_summary.csv`
- `results/paper2_revised3_multi_unknown_global.csv`
- `results/paper2_revised3_multi_unknown_by_scenario.csv`
- `results/paper2_revised3_ablation_summary.csv`
- `results/paper2_dt_runtime_audit.csv`
- `results/paper2_trigger_only_audit.csv`

## Suggested release notes

DT-RLR v1.1.0 adds the parameter-sensitivity, module-ablation, multiple-UNKNOWN, and computational-overhead experiments used in the revised manuscript targeted to Measurement. The core DT-RLR implementation and principal frozen configuration are unchanged.
