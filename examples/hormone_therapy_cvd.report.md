# SIGNET report: Hormone therapy and cardiovascular disease, evidence available in 1994

Input SHA-256: 5b41a98a867b06246741d60f8aad9e0e74f8fbe29c8c7e7a7d49741d9e5f3e72

## Target: cvd

| Rung | Capacity (upper, lower) | 95% SIGNET set (ratio scale) | k = 1 set | p benefit | p harm | Model test p |
| --- | --- | --- | --- | --- | --- | --- |
| rung1 | unbounded, bounded | 0.495 to inf | 0 to inf | 1.0000 | 1.0000 | 1.0 |
| rung2 | unbounded, bounded | 0.550 to inf | 0.495 to inf | 1.0000 | 1.0000 | 1.0 |
| rung3 | bounded, bounded | 0.550 to 0.858 | 0.495 to inf | 0.0003 | 1.0000 | 1.0 |

Default at rung3, tolerant to one failed assumption group or one broken evidence cell: 0 to inf
Beside it, rung3 without tolerance (conditional on every named assumption): 0.550 to 0.858
Capacity of the default: upper limit unbounded whatever the data; lower limit unbounded whatever the data.
One cell measures this target, so the default cannot support a directional claim here; its silence carries no information about the effect.

