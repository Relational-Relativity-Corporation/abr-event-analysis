"""
analyze_prediction_comparison.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Prediction test: which model is more correct?

A model is more correct when it makes predictions that are:
  1. More precise -- narrower bounds on what will be observed
  2. More complete -- accounts for more observable variance
  3. Traceable -- every predicted quantity traceable to observable through M

Method:
  Split 10,000 events: training (0-7999), test (8000-9999).

  QM prediction for each test event:
    invariant_mass = declared M_D0 (QM predicts the peak exactly)
    decay_time = declared tau (QM predicts the mean)
    kaon direction = (0,0,0) normalized -- QM predicts isotropic,
                     so the best single prediction is the mean = no direction
    kaon magnitude = p_star (fixed by 2-body kinematics -- QM predicts this)

  Operator prediction for each test event:
    Uses the relational field over the training set.
    For each test event at position N, the operator prediction is:
      x_predicted[N] = x[N-1] + mean(Delta field over training edges
                       incident to events near N-1)
    This is the simplest admissible operator prediction:
    the next observable state is the prior state plus the declared
    mean directed difference in its neighborhood.

  Comparison metric:
    Mean squared error (MSE) per component on the test set.
    Declared C projection:
      Preserves: prediction error per component.
      Discards: individual event values and ordering.

All targets declared before execution. No result adjusted post-hoc.
"""

import numpy as np
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from declaration.observable import OBSERVABLE_COMPONENTS, N_COMPONENTS
from declaration.relations import declare_relations, RelationProvenance
from declaration.targets import DECLARED_TARGETS
from operators.primary import operator_delta, operator_sigma, antisymmetric_term
from analysis.rank import im_delta_rank, SVD_TOLERANCE

EVENTS_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "data", "events", "d0_events.csv"
)
REPORT_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "results", "prediction_comparison_report.md"
)

# Declared split -- before execution
TRAIN_SIZE = 8000
TEST_SIZE  = 2000

# QM declared predictions (all PDG-traceable, declared before execution)
M_D0   = DECLARED_TARGETS["D0_MASS_MEV"]
TAU_D0 = DECLARED_TARGETS["D0_LIFETIME_PS"]
M_KAON = DECLARED_TARGETS["KAON_MASS_MEV"]
M_PION = DECLARED_TARGETS["PION_MASS_MEV"]

# Rest-frame kaon momentum magnitude -- fixed by 2-body kinematics
# This is a declared quantity computable from PDG masses
P_STAR = np.sqrt(
    (M_D0**2 - (M_KAON + M_PION)**2) *
    (M_D0**2 - (M_KAON - M_PION)**2)
) / (2 * M_D0)


def qm_prediction(n_test, D0_momentum_MeV=5000.0):
    """
    QM statistical prediction for each test event.

    Declared predictions (all traceable to PDG through M):
      invariant_mass: M_D0 (peak of Breit-Wigner -- QM predicts the mass)
      decay_time:     TAU_D0 (mean of exponential -- QM predicts the lifetime)
      decay_length:   derived from decay_time and declared beam momentum
      kaon momentum:  magnitude = P_STAR (fixed by kinematics)
                      direction = (0, 0, 1) normalized -- beam direction
                      QM predicts isotropic; best single prediction
                      is along beam (maximizes coverage, minimizes worst-case
                      error for isotropic distribution)
      vertex:         (0, 0, decay_length) -- along beam axis
      pion momentum:  equal and opposite to kaon along beam

    Declared C projection:
      Preserves: QM prediction per component.
      Discards: event-to-event variation (QM cannot predict individual events).
    """
    D0_gamma = np.sqrt((D0_momentum_MeV / M_D0)**2 + 1)
    D0_beta  = D0_momentum_MeV / (D0_gamma * M_D0)
    C_MM_PS  = 299.792458

    # Declared decay length from declared lifetime and beam
    declared_decay_length = TAU_D0 * D0_beta * C_MM_PS * D0_gamma

    # QM kaon momentum prediction: along beam (z), magnitude boosted
    E_kaon_rest = np.sqrt(M_KAON**2 + P_STAR**2)
    kaon_pz_lab = D0_gamma * (P_STAR + D0_beta * E_kaon_rest)
    E_pion_rest = np.sqrt(M_PION**2 + P_STAR**2)
    pion_pz_lab = D0_gamma * (-P_STAR + D0_beta * E_pion_rest)

    # Build prediction array -- same value for every test event
    pred = np.zeros((n_test, N_COMPONENTS))
    comp_idx = {c: i for i, c in enumerate(OBSERVABLE_COMPONENTS)}

    pred[:, comp_idx["vertex_x"]]      = 0.0
    pred[:, comp_idx["vertex_y"]]      = 0.0
    pred[:, comp_idx["vertex_z"]]      = declared_decay_length
    pred[:, comp_idx["kaon_px"]]       = 0.0
    pred[:, comp_idx["kaon_py"]]       = 0.0
    pred[:, comp_idx["kaon_pz"]]       = kaon_pz_lab
    pred[:, comp_idx["pion_px"]]       = 0.0
    pred[:, comp_idx["pion_py"]]       = 0.0
    pred[:, comp_idx["pion_pz"]]       = pion_pz_lab
    pred[:, comp_idx["decay_length"]]  = declared_decay_length
    pred[:, comp_idx["invariant_mass"]]= M_D0
    pred[:, comp_idx["decay_time"]]    = TAU_D0

    return pred


def operator_prediction(x_train, edges_train, x_test_prev):
    """
    Operator prediction for each test event.

    Method: for test event at index i (first test event = train[-1] + 1),
    the operator prediction is:

      x_pred[i] = x[i-1] + mean_delta

    where mean_delta is the mean directed difference across the training
    edge field -- the declared mean relational step over the training set.

    This is the simplest admissible operator prediction:
    the next observable state equals the prior state plus the mean
    directed difference observed in the training relational field.

    For components where the mean directed difference is near zero
    (e.g. invariant_mass, which is conserved), this reduces to x[i-1].
    For components where the mean directed difference is non-zero
    (e.g. decay_time, vertex_z), this extrapolates the trend.

    Declared C projection:
      Preserves: mean directed difference per component over training set.
      Discards: higher-order relational structure (Sigma output).
    """
    # Compute delta field over training set
    delta_field = operator_delta(x_train, edges_train)

    # Mean directed difference per component over training edges
    mean_delta = delta_field.mean(axis=0)  # shape (N_COMPONENTS,)

    # Std of directed differences -- measures how much variation exists
    std_delta = delta_field.std(axis=0)

    # Predict each test event as prior event + mean_delta
    n_test = len(x_test_prev)
    pred = np.zeros((n_test, N_COMPONENTS))
    for i, x_prev in enumerate(x_test_prev):
        pred[i] = x_prev + mean_delta

    return pred, mean_delta, std_delta


def operator_sigma_prediction(x_train, edges_train, x_test_prev,
                               test_event_ids, all_edges):
    """
    Enhanced operator prediction using Sigma antisymmetric term.

    For each test event, the prediction uses not just the mean delta
    but the local antisymmetric structure from Sigma:

      x_pred[i] = x[i-1] + delta_mean + sigma_correction[i-1]

    where sigma_correction captures the local asymmetric structure
    at the edge connecting event i-1 to event i.

    This is the Sigma-enhanced prediction -- it uses the full
    E_primary = Sigma(Delta(x)) output rather than just Delta mean.

    Declared C projection:
      Preserves: local antisymmetric structure per declared edge.
      Discards: global field structure beyond immediate adjacency.
    """
    delta_field = operator_delta(x_train, edges_train)
    sigma_field = operator_sigma(delta_field, edges_train, rho_base=0.2)

    # For each test event, use the Sigma output at the last training edge
    # as the correction term
    # Last training edge connects event TRAIN_SIZE-2 to TRAIN_SIZE-1
    last_sigma = sigma_field[-1]  # Sigma at the boundary edge

    n_test = len(x_test_prev)
    pred = np.zeros((n_test, N_COMPONENTS))
    for i, x_prev in enumerate(x_test_prev):
        pred[i] = x_prev + last_sigma

    return pred


def prediction_error(predicted, actual, label):
    """
    Compute prediction error per component.

    Declared C projection:
      Preserves: MSE per component, total MSE, relative improvement.
      Discards: individual event residuals.
    """
    residuals = predicted - actual
    mse_per_comp = (residuals**2).mean(axis=0)
    total_mse = mse_per_comp.mean()
    mae_per_comp = np.abs(residuals).mean(axis=0)

    print(f"\n  {label}:")
    print(f"    Total MSE:  {total_mse:.4f}")
    print(f"    Total MAE:  {mae_per_comp.mean():.4f}")
    print(f"\n    Per-component MAE:")
    for i, comp in enumerate(OBSERVABLE_COMPONENTS):
        print(f"      {comp:<20} {mae_per_comp[i]:.4f}")

    return mse_per_comp, total_mse, mae_per_comp


def run():
    print("\nanalyze_prediction_comparison.py")
    print("Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.")
    print("=" * 60)
    print()
    print("Prediction test: operator framework vs QM statistical model")
    print()
    print(f"Train set: events 0-{TRAIN_SIZE-1}  ({TRAIN_SIZE} events)")
    print(f"Test set:  events {TRAIN_SIZE}-{TRAIN_SIZE+TEST_SIZE-1}  ({TEST_SIZE} events)")
    print()
    print("Declared QM predictions:")
    print(f"  invariant_mass = {M_D0} MeV/c^2 (PDG mass)")
    print(f"  decay_time     = {TAU_D0} ps (PDG lifetime)")
    print(f"  kaon direction = along beam (z) -- isotropic mean")
    print(f"  kaon magnitude = {P_STAR:.4f} MeV/c (2-body kinematics)")
    print()

    # Load events
    if not os.path.exists(EVENTS_PATH):
        print(f"Event data not found: {EVENTS_PATH}")
        print("Run generate_events.py first.")
        sys.exit(1)

    df = pd.read_csv(EVENTS_PATH)
    n_total = len(df)
    print(f"Total events loaded: {n_total}")

    if n_total < TRAIN_SIZE + TEST_SIZE:
        print(f"Need {TRAIN_SIZE + TEST_SIZE} events. Got {n_total}.")
        print("Run: python generate_events.py --n 10000")
        sys.exit(1)

    # Build observable arrays
    X = np.array([[row[c] for c in OBSERVABLE_COMPONENTS]
                  for _, row in df.iterrows()])

    X_train = X[:TRAIN_SIZE]
    X_test  = X[TRAIN_SIZE:TRAIN_SIZE + TEST_SIZE]

    # Previous events for test prediction (each test event predicted from prior)
    X_test_prev = X[TRAIN_SIZE-1:TRAIN_SIZE + TEST_SIZE - 1]

    # Build training edge field
    x_train_dict = {i: X_train[i] for i in range(TRAIN_SIZE)}
    edges_train = declare_relations(
        list(x_train_dict.keys()),
        RelationProvenance.TEMPORAL,
        "declared production order -- training set"
    )

    # ── QM prediction ─────────────────────────────────────────────────────

    print("=" * 60)
    print("QM STATISTICAL MODEL PREDICTIONS")
    print()
    print("QM predicts distributions, not individual events.")
    print("Best single prediction = distribution mean for each component.")

    pred_qm = qm_prediction(TEST_SIZE)
    mse_qm, total_mse_qm, mae_qm = prediction_error(
        pred_qm, X_test, "QM prediction (distribution mean)")

    # ── Operator Delta prediction ─────────────────────────────────────────

    print()
    print("=" * 60)
    print("OPERATOR DELTA PREDICTION")
    print()
    print("Prediction: x[N] = x[N-1] + mean(Delta field over training set)")
    print("Uses directed differences across training events.")
    print("No distribution assumed. No statistical prior.")

    pred_op, mean_delta, std_delta = operator_prediction(
        x_train_dict, edges_train, X_test_prev)

    print(f"\n  Mean directed difference per component (training set):")
    for i, comp in enumerate(OBSERVABLE_COMPONENTS):
        print(f"    {comp:<20} mean={mean_delta[i]:+.4f}  std={std_delta[i]:.4f}")

    mse_op, total_mse_op, mae_op = prediction_error(
        pred_op, X_test, "Operator Delta prediction")

    # ── Operator Sigma prediction ─────────────────────────────────────────

    print()
    print("=" * 60)
    print("OPERATOR SIGMA PREDICTION (enhanced)")
    print()
    print("Prediction: x[N] = x[N-1] + Sigma output at boundary edge")
    print("Uses local antisymmetric structure from training field.")

    pred_sigma = operator_sigma_prediction(
        x_train_dict, edges_train, X_test_prev,
        list(range(TRAIN_SIZE, TRAIN_SIZE + TEST_SIZE)),
        edges_train
    )

    mse_sigma, total_mse_sigma, mae_sigma = prediction_error(
        pred_sigma, X_test, "Operator Sigma prediction")

    # ── Naive baseline ────────────────────────────────────────────────────

    print()
    print("=" * 60)
    print("NAIVE BASELINE (prior event -- no model)")
    print()
    print("Prediction: x[N] = x[N-1]  (no change)")
    print("This is the null model -- beats any model that adds noise.")

    mse_naive, total_mse_naive, mae_naive = prediction_error(
        X_test_prev, X_test, "Naive baseline (x[N] = x[N-1])")

    # ── Comparison ────────────────────────────────────────────────────────

    print()
    print("=" * 60)
    print("COMPARISON SUMMARY")
    print()
    print(f"  {'Model':<35} {'Total MSE':>12} {'Total MAE':>12}")
    print(f"  {'-'*35} {'-'*12} {'-'*12}")
    print(f"  {'Naive baseline (x[N]=x[N-1])':<35} "
          f"{total_mse_naive:>12.4f} {mae_naive.mean():>12.4f}")
    print(f"  {'QM (distribution mean)':<35} "
          f"{total_mse_qm:>12.4f} {mae_qm.mean():>12.4f}")
    print(f"  {'Operator Delta (mean relational step)':<35} "
          f"{total_mse_op:>12.4f} {mae_op.mean():>12.4f}")
    print(f"  {'Operator Sigma (local antisymmetric)':<35} "
          f"{total_mse_sigma:>12.4f} {mae_sigma.mean():>12.4f}")
    print()

    # Per-component comparison
    print(f"  Per-component MAE comparison:")
    print(f"  {'Component':<20} {'Naive':>10} {'QM':>10} "
          f"{'Op-Delta':>10} {'Op-Sigma':>10} {'Best':>12}")
    print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")

    results = []
    for i, comp in enumerate(OBSERVABLE_COMPONENTS):
        scores = {
            "Naive":    mae_naive[i],
            "QM":       mae_qm[i],
            "Op-Delta": mae_op[i],
            "Op-Sigma": mae_sigma[i],
        }
        best = min(scores, key=scores.get)
        results.append((comp, scores, best))
        print(f"  {comp:<20} {mae_naive[i]:>10.4f} {mae_qm[i]:>10.4f} "
              f"{mae_op[i]:>10.4f} {mae_sigma[i]:>10.4f} {best:>12}")

    print()

    # Count wins
    wins = {"Naive": 0, "QM": 0, "Op-Delta": 0, "Op-Sigma": 0}
    for _, scores, best in results:
        wins[best] += 1

    print(f"  Component wins:")
    for model, w in wins.items():
        print(f"    {model:<35} {w} / {N_COMPONENTS}")
    print()

    # Overall assessment
    op_wins_over_qm = sum(
        1 for _, scores, _ in results
        if min(scores["Op-Delta"], scores["Op-Sigma"]) < scores["QM"]
    )
    qm_wins_over_op = sum(
        1 for _, scores, _ in results
        if scores["QM"] < min(scores["Op-Delta"], scores["Op-Sigma"])
    )

    print(f"  Operator beats QM on: {op_wins_over_qm} / {N_COMPONENTS} components")
    print(f"  QM beats operator on: {qm_wins_over_op} / {N_COMPONENTS} components")
    print()

    if op_wins_over_qm > qm_wins_over_op:
        finding = "OPERATOR FRAMEWORK MORE PREDICTIVE"
        detail = ("The operator framework produces lower prediction error "
                  "than the QM statistical model on more observable components. "
                  "More information produces a more complete model.")
    elif op_wins_over_qm == qm_wins_over_op:
        finding = "EQUIVALENT PREDICTIVE POWER"
        detail = ("Operator framework and QM statistical model produce "
                  "equivalent prediction error. The additional rank dimensions "
                  "detected by the operators do not improve individual event "
                  "prediction for these components.")
    else:
        finding = "QM MORE PREDICTIVE ON THESE COMPONENTS"
        detail = ("QM statistical model produces lower prediction error "
                  "on more components. The operator prediction method "
                  "requires refinement.")

    print(f"  Finding: {finding}")
    print(f"  {detail}")
    print()

    # Key insight regardless of outcome
    print("=" * 60)
    print("KEY FINDING (independent of prediction comparison):")
    print()
    print("  The rank analysis establishes:")
    rank_train = im_delta_rank(operator_delta(x_train_dict, edges_train))
    print(f"  rank(Im Delta) over {TRAIN_SIZE} training events: {rank_train[0]}")
    print(f"  PDG summary rank: 3")
    print(f"  Information gap: {rank_train[0] - 3} dimensions")
    print()
    print("  Statistical reduction discards real relational structure.")
    print("  Whether that structure improves prediction depends on")
    print("  the declared prediction method -- which is an open question.")
    print("  The simplest operator prediction (mean delta) may not be")
    print("  the most powerful use of the relational field.")

    # Write report
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("# Prediction Comparison Report\n\n")
        f.write("**Metatron Dynamics, Inc.** V6. "
                "Bounded over D. No claim beyond D.\n\n")
        f.write("## Declaration\n\n")
        f.write(f"Train: {TRAIN_SIZE} events. Test: {TEST_SIZE} events.\n\n")
        f.write(f"QM predicted mass: {M_D0} MeV/c^2\n\n")
        f.write(f"QM predicted lifetime: {TAU_D0} ps\n\n")
        f.write("## Results\n\n")
        f.write(f"| Model | Total MSE | Total MAE |\n")
        f.write(f"|---|---|---|\n")
        f.write(f"| Naive | {total_mse_naive:.4f} | {mae_naive.mean():.4f} |\n")
        f.write(f"| QM | {total_mse_qm:.4f} | {mae_qm.mean():.4f} |\n")
        f.write(f"| Op-Delta | {total_mse_op:.4f} | {mae_op.mean():.4f} |\n")
        f.write(f"| Op-Sigma | {total_mse_sigma:.4f} | {mae_sigma.mean():.4f} |\n\n")
        f.write(f"Operator beats QM: {op_wins_over_qm}/{N_COMPONENTS} components\n\n")
        f.write(f"QM beats operator: {qm_wins_over_op}/{N_COMPONENTS} components\n\n")
        f.write(f"Finding: {finding}\n\n")
        f.write(f"rank(Im Delta) training set: {rank_train[0]}\n\n")
        f.write("## Declared projection\n\n")
        f.write("Preserves: MSE and MAE per component on held-out test set.\n\n")
        f.write("Discards: individual event residuals and ordering.\n\n")

    print()
    print(f"Report written: {REPORT_PATH}")
    print("\nBounded over D. No claim beyond D.")


if __name__ == "__main__":
    run()
