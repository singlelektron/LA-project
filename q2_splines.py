"""Question 2: cubic spline interpolation for atmospheric CO2 data.

The splines are built only from the seven data points printed in the PDF.
The NASA values for 1990 and 2010 are used afterward as an accuracy check.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib")))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache")))
Path(os.environ["MPLCONFIGDIR"]).mkdir(exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


YEARS = np.array([1850, 1875, 1900, 1925, 1950, 1975, 2000], dtype=float)
CO2_PPM = np.array([285.2, 288.6, 295.7, 305.3, 311.3, 331.36, 369.64], dtype=float)
TARGET_YEARS = np.array([1990.0, 2010.0], dtype=float)

# These values are used only after interpolation, for accuracy comparison.
ACTUAL_CO2_PPM = {
    1990: 354.29,
    2010: 389.21,
}

FIGURE_DIR = Path("figures")
POST_PATH = Path("post.md")


@dataclass(frozen=True)
class SplineResult:
    name: str
    slug: str
    moments: np.ndarray
    estimates: dict[int, float]


def build_rhs(years: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Build the right side for equally spaced cubic spline moment equations."""
    h = years[1] - years[0]
    second_differences = values[2:] - 2.0 * values[1:-1] + values[:-2]
    return 6.0 * second_differences / h**2


def interior_matrix(size: int) -> np.ndarray:
    """Matrix with diagonal 4 and off-diagonal 1 for interior moments."""
    matrix = np.zeros((size, size), dtype=float)
    np.fill_diagonal(matrix, 4.0)
    np.fill_diagonal(matrix[1:], 1.0)
    np.fill_diagonal(matrix[:, 1:], 1.0)
    return matrix


def natural_moments(years: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Second derivative values for the natural cubic spline."""
    rhs = build_rhs(years, values)
    matrix = interior_matrix(len(rhs))
    interior = np.linalg.solve(matrix, rhs)
    return np.concatenate(([0.0], interior, [0.0]))


def parabolic_runout_moments(years: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Second derivative values for the parabolic runout spline."""
    rhs = build_rhs(years, values)
    matrix = interior_matrix(len(rhs))
    matrix[0, 0] = 5.0
    matrix[-1, -1] = 5.0
    interior = np.linalg.solve(matrix, rhs)
    return np.concatenate(([interior[0]], interior, [interior[-1]]))


def cubic_runout_moments(years: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Second derivative values for the cubic runout spline."""
    rhs = build_rhs(years, values)
    matrix = interior_matrix(len(rhs))
    matrix[0, 0] = 6.0
    matrix[0, 1] = 0.0
    matrix[-1, -1] = 6.0
    matrix[-1, -2] = 0.0
    interior = np.linalg.solve(matrix, rhs)

    left_endpoint = 2.0 * interior[0] - interior[1]
    right_endpoint = 2.0 * interior[-1] - interior[-2]
    return np.concatenate(([left_endpoint], interior, [right_endpoint]))


def evaluate_spline(
    query_years: np.ndarray | list[float] | float,
    years: np.ndarray,
    values: np.ndarray,
    moments: np.ndarray,
) -> np.ndarray:
    """Evaluate a cubic spline, using the endpoint interval for extrapolation."""
    scalar_input = np.isscalar(query_years)
    query = np.atleast_1d(np.asarray(query_years, dtype=float))
    output = np.empty_like(query, dtype=float)

    for index, x in enumerate(query):
        if x <= years[0]:
            interval = 0
        elif x >= years[-1]:
            interval = len(years) - 2
        else:
            interval = int(np.searchsorted(years, x) - 1)

        x_left = years[interval]
        x_right = years[interval + 1]
        y_left = values[interval]
        y_right = values[interval + 1]
        m_left = moments[interval]
        m_right = moments[interval + 1]
        h = x_right - x_left

        a = (x_right - x) / h
        b = (x - x_left) / h
        output[index] = (
            a * y_left
            + b * y_right
            + ((a**3 - a) * m_left + (b**3 - b) * m_right) * h**2 / 6.0
        )

    if scalar_input:
        return np.array(output[0])
    return output


def compute_results() -> list[SplineResult]:
    """Compute all three spline variants and their target-year estimates."""
    methods = [
        ("Natural cubic spline", "natural_cubic_spline", natural_moments),
        ("Parabolic runout spline", "parabolic_runout_spline", parabolic_runout_moments),
        ("Cubic runout spline", "cubic_runout_spline", cubic_runout_moments),
    ]

    results = []
    for name, slug, moment_builder in methods:
        moments = moment_builder(YEARS, CO2_PPM)
        estimates = {
            int(year): float(evaluate_spline(year, YEARS, CO2_PPM, moments))
            for year in TARGET_YEARS
        }
        results.append(SplineResult(name=name, slug=slug, moments=moments, estimates=estimates))
    return results


def verify_interpolation(results: list[SplineResult]) -> None:
    """Ensure each spline reproduces all original data points."""
    for result in results:
        fitted = evaluate_spline(YEARS, YEARS, CO2_PPM, result.moments)
        if not np.allclose(fitted, CO2_PPM, atol=1e-9):
            raise AssertionError(f"{result.name} does not pass through the input data")


def plot_single_spline(result: SplineResult) -> None:
    """Save one plot for one spline method."""
    x_plot = np.linspace(YEARS[0], 2010.0, 500)
    y_plot = evaluate_spline(x_plot, YEARS, CO2_PPM, result.moments)

    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.plot(x_plot, y_plot, color="#1f77b4", linewidth=2.2, label=result.name)
    ax.scatter(YEARS, CO2_PPM, color="#111111", s=35, zorder=3, label="PDF data")
    ax.scatter(
        TARGET_YEARS,
        [result.estimates[int(year)] for year in TARGET_YEARS],
        color="#d62728",
        marker="x",
        s=70,
        zorder=4,
        label="Spline estimates",
    )
    ax.axvline(2000, color="#888888", linestyle="--", linewidth=1.0, alpha=0.75)
    ax.text(2001.5, min(y_plot) + 4, "extrapolation", color="#555555", fontsize=9)
    ax.set_title(result.name)
    ax.set_xlabel("Year")
    ax.set_ylabel("CO2 concentration (ppm)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / f"{result.slug}.png", dpi=200)
    plt.close(fig)


def plot_comparison(results: list[SplineResult]) -> None:
    """Save one plot comparing all spline methods."""
    x_plot = np.linspace(YEARS[0], 2010.0, 600)
    colors = {
        "Natural cubic spline": "#1f77b4",
        "Parabolic runout spline": "#2ca02c",
        "Cubic runout spline": "#d62728",
    }

    fig, ax = plt.subplots(figsize=(9, 5.4))
    for result in results:
        ax.plot(
            x_plot,
            evaluate_spline(x_plot, YEARS, CO2_PPM, result.moments),
            linewidth=2.0,
            color=colors[result.name],
            label=result.name,
        )

    ax.scatter(YEARS, CO2_PPM, color="#111111", s=35, zorder=3, label="PDF data")
    ax.scatter(
        TARGET_YEARS,
        [ACTUAL_CO2_PPM[int(year)] for year in TARGET_YEARS],
        color="#ff7f0e",
        marker="D",
        s=50,
        zorder=4,
        label="Actual NASA values",
    )
    ax.axvline(2000, color="#888888", linestyle="--", linewidth=1.0, alpha=0.75)
    ax.text(2001.5, 287, "2010 is extrapolated", color="#555555", fontsize=9)
    ax.set_title("Comparison of the three cubic spline interpolations")
    ax.set_xlabel("Year")
    ax.set_ylabel("CO2 concentration (ppm)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "spline_comparison.png", dpi=200)
    plt.close(fig)


def result_rows(results: list[SplineResult]) -> list[tuple[str, int, float, float, float, float]]:
    """Return rows containing estimates and errors."""
    rows = []
    for result in results:
        for year in sorted(ACTUAL_CO2_PPM):
            estimate = result.estimates[year]
            actual = ACTUAL_CO2_PPM[year]
            error = estimate - actual
            rows.append((result.name, year, estimate, actual, error, abs(error)))
    return rows


def print_results(results: list[SplineResult]) -> None:
    """Print a compact numerical result table."""
    print("Question 2 spline estimates")
    print("Input data: only the seven PDF table points")
    print("Actual values: used afterward for accuracy checking")
    print()
    print(f"{'Method':28s} {'Year':>6s} {'Estimate':>12s} {'Actual':>10s} {'Error':>11s} {'Abs error':>11s}")
    print("-" * 84)
    for method, year, estimate, actual, error, abs_error in result_rows(results):
        print(
            f"{method:28s} {year:6d} {estimate:12.4f} {actual:10.4f} "
            f"{error:11.4f} {abs_error:11.4f}"
        )


def write_post(results: list[SplineResult]) -> None:
    """Write the after-work explanation and usage notes."""
    rows = result_rows(results)
    table_lines = [
        "| Method | Year | Estimate (ppm) | Actual NASA value (ppm) | Error (ppm) | Absolute error (ppm) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for method, year, estimate, actual, error, abs_error in rows:
        table_lines.append(
            f"| {method} | {year} | {estimate:.4f} | {actual:.4f} | "
            f"{error:.4f} | {abs_error:.4f} |"
        )

    by_year = {}
    for method, year, estimate, actual, error, abs_error in rows:
        by_year.setdefault(year, []).append((method, abs_error))

    best_lines = []
    for year in sorted(by_year):
        best_method, best_error = min(by_year[year], key=lambda item: item[1])
        best_lines.append(
            f"- For {year}, the smallest absolute error is from the {best_method} "
            f"with absolute error {best_error:.4f} ppm."
        )

    POST_PATH.write_text(
        "\n".join(
            [
                "# Question 2 Work Summary",
                "",
                "## What Was Done",
                "",
                "I implemented the three cubic spline interpolation methods requested in Question 2:",
                "",
                "- Natural cubic spline",
                "- Parabolic runout spline",
                "- Cubic runout spline",
                "",
                "The splines are built only from the seven CO2 data points printed in the PDF. "
                "The NASA values are used afterward only to check accuracy.",
                "",
                "## Files Created",
                "",
                "- `q2_splines.py`: Python implementation for all three spline methods.",
                "- `figures/natural_cubic_spline.png`: Plot of the natural cubic spline.",
                "- `figures/parabolic_runout_spline.png`: Plot of the parabolic runout spline.",
                "- `figures/cubic_runout_spline.png`: Plot of the cubic runout spline.",
                "- `figures/spline_comparison.png`: One plot comparing all three splines.",
                "- `post.md`: This summary and usage note.",
                "",
                "## How To Run",
                "",
                "From the project directory, run:",
                "",
                "```bash",
                ".venv/bin/python q2_splines.py",
                "```",
                "",
                "The script prints the numerical result table, regenerates the plots, and rewrites `post.md`.",
                "",
                "## Numerical Results",
                "",
                *table_lines,
                "",
                "## Comments On Accuracy",
                "",
                *best_lines,
                "",
                "The 1990 value is an interpolation because it lies between 1975 and 2000. "
                "The 2010 value is an extrapolation because it lies beyond the last PDF data point, 2000, "
                "so it should be treated as less reliable.",
                "",
                "For this data set, the natural spline is closer for 1990, while the cubic runout spline "
                "is closer for 2010. The parabolic and cubic runout splines extrapolate more realistically "
                "past 2000 than the natural spline because they do not force the final endpoint curvature to zero.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    results = compute_results()
    verify_interpolation(results)
    for result in results:
        plot_single_spline(result)
    plot_comparison(results)
    print_results(results)
    write_post(results)


if __name__ == "__main__":
    main()
