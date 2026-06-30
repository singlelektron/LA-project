# Question 2 Work Plan

## What Question 2 Asks

Question 2 asks us to use cubic spline interpolation on the atmospheric CO2 data printed in the project PDF. The spline must be built from the table in the PDF, then used to estimate CO2 concentration for 1990 and 2010.

After the spline estimates are computed, the NASA data link from the PDF is used only as a check of accuracy. The outside data is not used to construct any spline.

## Data Used For Interpolation

The interpolation data from the PDF is:

| Year | CO2 (ppm) |
| ---: | ---: |
| 1850 | 285.20 |
| 1875 | 288.60 |
| 1900 | 295.70 |
| 1925 | 305.30 |
| 1950 | 311.30 |
| 1975 | 331.36 |
| 2000 | 369.64 |

The years are equally spaced, so

```text
h = 25 years
```

## Mathematical Setup

Let the data points be `(x_i, y_i)`, where `x_i` is the year and `y_i` is the CO2 value. The cubic spline is written using the second derivative values

```text
M_i = S''(x_i)
```

For each interval `[x_i, x_{i+1}]`, the spline value is computed by

```text
S_i(x) =
  M_i (x_{i+1} - x)^3 / (6h)
  + M_{i+1} (x - x_i)^3 / (6h)
  + (y_i - M_i h^2 / 6) (x_{i+1} - x) / h
  + (y_{i+1} - M_{i+1} h^2 / 6) (x - x_i) / h
```

The unknown interior second derivatives are found by solving a tridiagonal linear system. Since the data is equally spaced, the right side is

```text
6 / h^2 * (y_{i+1} - 2y_i + y_{i-1})
```

for the interior data points.

## Boundary Conditions To Implement

### 1. Natural Cubic Spline

The natural spline assumes zero curvature at both endpoints:

```text
M_0 = 0
M_n = 0
```

The interior system has diagonal entries `4` and off-diagonal entries `1`.

### 2. Parabolic Runout Spline

The parabolic runout spline assumes the first two second derivatives are equal, and the last two second derivatives are equal:

```text
M_0 = M_1
M_n = M_{n-1}
```

This changes the first and last diagonal entries of the interior system from `4` to `5`.

### 3. Cubic Runout Spline

The cubic runout spline uses endpoint third-derivative runout conditions. For equally spaced data, this gives

```text
M_0 = 2M_1 - M_2
M_n = 2M_{n-1} - M_{n-2}
```

The first and last rows of the interior system have diagonal entry `6` and no adjacent endpoint coefficient.

## Planned Output Files

The Python program will create:

```text
q2_splines.py
figures/natural_cubic_spline.png
figures/parabolic_runout_spline.png
figures/cubic_runout_spline.png
figures/spline_comparison.png
post.md
```

The script will print a table of estimates and errors. It will also save the same numerical results in `post.md`.

## Comparison Rule

The comparison values for 1990 and 2010 come from the NASA link in the PDF. They are used only after all spline estimates have already been calculated from the seven PDF data points.
