# Gold implied vol: GVZ vs MOVE vs VIX

Tests whether gold implied volatility (GVZ) co-moves more with macro/rates risk (MOVE) or equity volatility (VIX) during a GLD uptrend, combining descriptive evidence and out-of-sample volatility model comparison (QLIKE + MCS).

## Question
During Jun 2025 – Jan 19, 2026, does GVZ track MOVE more than VIX while GLD is trending up?

## Data
Pulled from Yahoo Finance via `quantmod`:
- GVZ (CBOE Gold ETF Volatility Index)
- VIX (CBOE Volatility Index)
- MOVE (ICE BofA MOVE Index)
- OVX (CBOE Crude Oil ETF Volatility Index)
- GLD (SPDR Gold Shares)

## Methods
- Rolling correlations: 30D and 60D (GVZ vs MOVE, GVZ vs VIX)
- Lagged regression on levels: `GVZ ~ MOVE(t-1) + VIX(t-1)`
- Volatility forecasting on GLD returns:
  - Models: sGARCH, eGARCH, apARCH (+ optional EGARCH-X with external regressors)
  - Metric: QLIKE
  - Model Confidence Set (MCS) out-of-sample

## Repository structure
- `R/` core functions (pipelines, plots, utilities)
- `scripts/` runnable entry point
- `reports/` report helpers (if used)

## How to run
From the repository root:

```bash
Rscript scripts/run_all.R
