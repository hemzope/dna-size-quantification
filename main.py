import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --------------- CREATE FOLDERS ---------------

RESULTS_DIR = "results"

LADDER_DIR = os.path.join(RESULTS_DIR, "ladder")
EVIDENCE_DIR = os.path.join(RESULTS_DIR, "evidence")
SUSPECTS_DIR = os.path.join(RESULTS_DIR, "suspects")
COMPARISON_DIR = os.path.join(RESULTS_DIR, "comparison")

for directory in [
    LADDER_DIR,
    EVIDENCE_DIR,
    SUSPECTS_DIR,
    COMPARISON_DIR
]:
    os.makedirs(directory, exist_ok=True)

# --------------- INPUT DATA ---------------

# ladder_bp = np.array([
#     3000, 1500, 1000, 900,
#     800, 700, 600, 500,
#     400, 300, 200
# ])

# ladder_distance = np.array([
#     262, 397, 518, 552, 590, 624, 682, 
#     740, 804, 880, 965
# ])

# evidence_distance = np.array([
#     365, 465, 585, 828
# ])

# suspect1_distance = np.array([
#     272, 373, 430, 730, 833
# ])

# suspect2_distance = np.array([
#     335, 445, 573, 822
# ])

# suspect3_distance = np.array([
#     340, 445, 570, 610, 827
# ])

with open("dna_data.json", "r") as f:
    data = json.load(f)

ladder_bp = np.array(
    data["ladder"]["bp"]
)

ladder_distance = np.array(
    data["ladder"]["distance"]
)

evidence_distance = np.array(
    data["evidence"]["distance"]
)

suspects = {}

for suspect_name, suspect_data in data["suspects"].items():
    suspects[suspect_name] = np.array(
        suspect_data["distance"]
    )

# ============================================================
# VALIDATION
# ============================================================

if len(ladder_bp) != len(ladder_distance):
    raise ValueError(
        "ladder_bp and ladder_distance must have the same number of values."
    )

# ============================================================
# STANDARD CURVE
# ============================================================

ladder_log_bp = np.log10(ladder_bp)

slope, intercept = np.polyfit(
    ladder_distance,
    ladder_log_bp,
    1
)

predicted_ladder = slope * ladder_distance + intercept

ss_res = np.sum(
    (ladder_log_bp - predicted_ladder) ** 2
)

ss_tot = np.sum(
    (ladder_log_bp - np.mean(ladder_log_bp)) ** 2
)

r2 = 1 - (ss_res / ss_tot)

print("\n===== STANDARD CURVE =====")
print(f"Slope      : {slope:.6f}")
print(f"Intercept  : {intercept:.6f}")
print(f"R²         : {r2:.6f}")

# ============================================================
# SAVE LADDER CSV
# ============================================================

ladder_df = pd.DataFrame({
    "Distance": ladder_distance,
    "Base_Pairs": ladder_bp,
    "Log10_BP": ladder_log_bp
})

ladder_df.to_csv(
    os.path.join(LADDER_DIR, "ladder_data.csv"),
    index=False
)

# ============================================================
# FUNCTIONS
# ============================================================

def estimate_bp(distances):
    """
    Convert migration distances into estimated bp values.
    """

    log_bp = slope * distances + intercept
    bp = 10 ** log_bp

    return log_bp, bp


def save_sample_csv(name, distances, folder):
    """
    Save estimated DNA sizes.
    """

    log_bp, bp = estimate_bp(distances)

    df = pd.DataFrame({
        "Distance": distances,
        "Estimated_Log10_BP": log_bp,
        "Estimated_BP": bp
    })

    filepath = os.path.join(
        folder,
        f"{name}_results.csv"
    )

    df.to_csv(filepath, index=False)

    return log_bp, bp


def create_plot(name, distances, folder):
    """
    Create DNA quantification plot.
    """

    unknown_log_bp, _ = estimate_bp(distances)

    x_fit = np.linspace(
        ladder_distance.min(),
        ladder_distance.max(),
        300
    )

    y_fit = slope * x_fit + intercept

    plt.figure(figsize=(8, 6))

    # Ladder points
    plt.scatter(
        ladder_distance,
        ladder_log_bp,
        label="DNA Ladder",
        color='black'
    )

    # Trendline
    plt.plot(
        x_fit,
        y_fit,
        linewidth=2,
        label="Trendline",
        color='red'
    )

    # Unknown sample
    plt.scatter(
        distances,
        unknown_log_bp,
        marker='x',
        s=120,
        label=name.capitalize()
    )

    equation = (
        f"log10(bp) = {slope:.6f}x + {intercept:.6f}\n"
        f"R² = {r2:.5f}"
    )

    plt.text(
        0.86,
        0.7,
        equation,
        transform=plt.gca().transAxes,
        horizontalalignment="right",
        verticalalignment="top",
        color="black",
        bbox=dict(
            boxstyle="round",
            facecolor="lightgray",
            edgecolor="black",
            alpha=0.2
        )
    )

    plt.xlabel("Migration Distance")
    plt.ylabel("log10(Base Pairs)")
    plt.title(f"DNA Size Quantification - {name.capitalize()}")

    plt.grid(
        True,
        linestyle='--',
        alpha=0.7
    )
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        os.path.join(folder, f"{name}_plot.png"),
        dpi=300
    )

    plt.close()


def matching_error(evidence_log_bp, suspect_log_bp):
    """
    Mean absolute error in log-space.
    """

    n = min(
        len(evidence_log_bp),
        len(suspect_log_bp)
    )

    if n == 0:
        return np.inf

    evidence_log_bp = evidence_log_bp[:n]
    suspect_log_bp = suspect_log_bp[:n]

    return np.mean(
        np.abs(
            evidence_log_bp - suspect_log_bp
        )
    )

# ============================================================
# SAVE STANDARD CURVE PLOT
# ============================================================

x_fit = np.linspace(
    ladder_distance.min(),
    ladder_distance.max(),
    300
)

y_fit = slope * x_fit + intercept

plt.figure(figsize=(8, 6))

plt.scatter(
    ladder_distance,
    ladder_log_bp,
    label="DNA Ladder",
    color='black'
)

plt.plot(
    x_fit,
    y_fit,
    linewidth=2,
    label="Trendline",
    color='red'
)

equation = (
    f"R² = {r2:.5f}\n"
    f"log10(bp) = {slope:.6f}x + {intercept:.6f}"
)

plt.text(
    0.86,
    0.7,
    equation,
    transform=plt.gca().transAxes,
    horizontalalignment="right",
    verticalalignment="top",
    color="black",
    bbox=dict(
        boxstyle="round",
        facecolor="lightgray",
        edgecolor="black",
        alpha=0.2
    )
)

plt.xlabel("Migration Distance")
plt.ylabel("log10(Base Pairs)")
plt.title("DNA Ladder Standard Curve")

plt.grid(
    True,
    linestyle='--',
    alpha=0.7
)
plt.legend()
plt.tight_layout()

plt.savefig(
    os.path.join(
        LADDER_DIR,
        "standard_curve.png"
    ),
    dpi=300
)

plt.close()

# ============================================================
# PROCESS EVIDENCE
# ============================================================

evidence_log_bp, evidence_bp = save_sample_csv(
    "evidence",
    evidence_distance,
    EVIDENCE_DIR
)

create_plot(
    "evidence",
    evidence_distance,
    EVIDENCE_DIR
)

# ============================================================
# PROCESS SUSPECTS
# ============================================================

scores = []

for name, distances in suspects.items():

    suspect_log_bp, suspect_bp = save_sample_csv(
        name,
        distances,
        SUSPECTS_DIR
    )

    create_plot(
        name,
        distances,
        SUSPECTS_DIR
    )

    score = matching_error(
        evidence_log_bp,
        suspect_log_bp
    )

    scores.append(
        [name, score]
    )

# ============================================================
# RANK SUSPECTS
# ============================================================

ranking_df = pd.DataFrame(
    scores,
    columns=[
        "Suspect",
        "Mean_Absolute_Log10_Error"
    ]
)

ranking_df = ranking_df.sort_values(
    by="Mean_Absolute_Log10_Error"
)

ranking_df.to_csv(
    os.path.join(
        COMPARISON_DIR,
        "suspect_ranking.csv"
    ),
    index=False
)

best_match = ranking_df.iloc[0]

# ============================================================
# MATCH REPORT
# ============================================================

report_df = pd.DataFrame({
    "Slope": [slope],
    "Intercept": [intercept],
    "R2": [r2],
    "Number_of_Suspects": len(suspects),
    "Best_Match": [best_match["Suspect"]],
    "Matching_Error": [
        best_match["Mean_Absolute_Log10_Error"]
    ]
})

report_df.to_csv(
    os.path.join(
        COMPARISON_DIR,
        "match_report.csv"
    ),
    index=False
)

# ============================================================
# PRINT SUMMARY
# ============================================================

print("\n===== SUSPECT RANKING =====")

for _, row in ranking_df.iterrows():
    print(
        f"{row['Suspect']}: "
        f"{row['Mean_Absolute_Log10_Error']:.6f}"
    )

print(
    f"\nBEST MATCH: {best_match['Suspect']}"
)

print("\nAnalysis complete.")

print(f"\nResults saved in:")
print(f"  {RESULTS_DIR}/")