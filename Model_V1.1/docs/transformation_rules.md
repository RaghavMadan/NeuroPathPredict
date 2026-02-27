# Covariate Transformation Rules

This document describes the rules used to suggest transformations for covariates in the NeuroPathPredict pipeline.

## Overview
The pipeline analyzes each covariate's distribution characteristics and suggests appropriate transformations based on:
- Skewness
- Kurtosis
- Range (presence of negative values)
- Variance
- Unique value count

## Transformation Rules

### 1. Constant/Near-Constant Values
- **Rule**: Standard deviation < 1e-10
- **Suggestion**: "constant"
- **Rationale**: These variables have essentially no variation and might be removed from analysis

### 2. Binary/Categorical Data
- **Rule**: Number of unique values ≤ 2 (excluding zeros)
- **Suggestion**: "binary"
- **Rationale**: These variables are already in a categorical format and shouldn't be transformed

### 3. Highly Skewed Positive Data
- **Rule**: Minimum value ≥ 0 AND |skewness| > 3
- **Suggestion**: "log1p"
- **Rationale**: Log transformation (log(x + 1)) is effective for highly right-skewed positive data
- **Note**: Using log1p instead of log to handle zeros

### 4. Moderately Skewed Positive Data
- **Rule**: Minimum value ≥ 0 AND |skewness| > 1.5
- **Suggestion**: "sqrt"
- **Rationale**: Square root transformation is less aggressive than log and suitable for moderate skewness

### 5. Skewed Data with Negative Values
- **Rule**: Minimum value < 0
- Two sub-cases:
  a. **Extremely Skewed** (|skewness| > 3)
     - **Suggestion**: "cbrt" (cube root)
     - **Rationale**: Cube root is effective for highly skewed data with negative values
  b. **Moderately Skewed** (|skewness| > 1.5)
     - If heavy tails (|kurtosis| > 3):
       - **Suggestion**: "cbrt" (cube root)
       - **Rationale**: Cube root handles negative values while reducing skewness and heavy tails
     - If normal tails:
       - **Suggestion**: "boxcox"
       - **Rationale**: Box-Cox transformation can handle various shapes after data shifting

### 6. Symmetric Heavy-Tailed Data
- **Rule**: |skewness| ≤ 1.5 AND |kurtosis| > 3
- **Suggestion**: "scale"
- **Rationale**: Data is symmetric but has heavy tails; standardization might be sufficient

### 7. Well-Behaved Data
- **Rule**: None of the above conditions met
- **Suggestion**: "none"
- **Rationale**: Data is reasonably well-behaved and doesn't require transformation

## Implementation Notes

1. Transformations are applied in order of specificity:
   - Check for constants/near-constants first
   - Then check for binary/categorical data
   - Then check for skewness patterns
   - Finally check for heavy tails

2. Thresholds used:
   - Severe skewness: |skewness| > 3
   - Moderate skewness: |skewness| > 1.5
   - Heavy tails: |kurtosis| > 3
   - Near-constant: sd < 1e-10

3. After transformation, it's recommended to:
   - Verify the transformation improved the distribution
   - Check for any introduced artifacts
   - Consider the interpretability of the transformed variable

## References
- Skewness thresholds based on Bulmer (1979)
- Kurtosis thresholds based on DeCarlo (1997)
- Transformation strategies adapted from Box-Cox (1964) and Tukey (1977) 