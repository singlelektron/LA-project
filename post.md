# Question 2 Work Summary

## What Was Done

I implemented the three cubic spline interpolation methods requested in Question 2:

- Natural cubic spline
- Parabolic runout spline
- Cubic runout spline

The splines are built only from the seven CO2 data points printed in the PDF. The NASA values are used afterward only to check accuracy.

## Files Created

- `q2_splines.py`: Python implementation for all three spline methods.
- `figures/natural_cubic_spline.png`: Plot of the natural cubic spline.
- `figures/parabolic_runout_spline.png`: Plot of the parabolic runout spline.
- `figures/cubic_runout_spline.png`: Plot of the cubic runout spline.
- `figures/spline_comparison.png`: One plot comparing all three splines.
- `post.md`: This summary and usage note.

## How To Run

From the project directory, run:

```bash
.venv/bin/python q2_splines.py
```

The script prints the numerical result table, regenerates the plots, and rewrites `post.md`.

## Numerical Results

| Method | Year | Estimate (ppm) | Actual NASA value (ppm) | Error (ppm) | Absolute error (ppm) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Natural cubic spline | 1990 | 353.0528 | 354.2900 | -1.2372 | 1.2372 |
| Natural cubic spline | 2010 | 386.2272 | 389.2100 | -2.9828 | 2.9828 |
| Parabolic runout spline | 1990 | 352.1733 | 354.2900 | -2.1167 | 2.1167 |
| Parabolic runout spline | 2010 | 389.9797 | 389.2100 | 0.7697 | 0.7697 |
| Cubic runout spline | 1990 | 352.2210 | 354.2900 | -2.0690 | 2.0690 |
| Cubic runout spline | 2010 | 389.7757 | 389.2100 | 0.5657 | 0.5657 |

## Comments On Accuracy

- For 1990, the smallest absolute error is from the Natural cubic spline with absolute error 1.2372 ppm.
- For 2010, the smallest absolute error is from the Cubic runout spline with absolute error 0.5657 ppm.

The 1990 value is an interpolation because it lies between 1975 and 2000. The 2010 value is an extrapolation because it lies beyond the last PDF data point, 2000, so it should be treated as less reliable.

For this data set, the natural spline is closer for 1990, while the cubic runout spline is closer for 2010. The parabolic and cubic runout splines extrapolate more realistically past 2000 than the natural spline because they do not force the final endpoint curvature to zero.
