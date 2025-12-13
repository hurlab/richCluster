## Resubmission - richCluster 1.0.2

This is a resubmission addressing all feedback from the CRAN review.

### Changes Made

1. **Documentation - Missing \value tag**:
   - Added comprehensive `\value` documentation to `runRichCluster()` function in R/cluster.R
   - The return value is now fully documented as a list containing `distance_matrix`, `all_clusters`, and `linkage_tree` components with detailed descriptions of each element

2. **Examples - \dontrun{} vs \donttest{}**:
   - Replaced `\dontrun{}` with `\donttest{}` in all example sections for functions: `term_bar()`, `term_dot()`, `cluster_bar()`, and `cluster_dot()` in R/bars.R and R/dots.R
   - Added code to load example data (`cluster_result.rds` from package extdata) in each example
   - These examples are now fully executable using the package's built-in example data, making `\donttest{}` the appropriate wrapper

3. **Console Output - Unsuppressible messages**:
   - Removed `print(title)` statement at line 196 in R/hmaps.R (term_hmap function)
   - This ensures all console output can be properly suppressed by users

4. **Graphics Parameters - par() settings**:
   - Added `on.exit(graphics::par(oldpar))` in R/plot_network_graph.R at line 94
   - User's graphics parameters are now properly saved before modification and restored upon function exit, even if the function breaks
   - The `on.exit()` call is placed immediately after saving the original settings as recommended

All issues identified in the previous review have been fully addressed.

## Test environments

* Local: Windows 11, R 4.4.2
* GitHub Actions (via R-CMD-check):
  - Windows (latest), R release
  - macOS (latest), R release
  - Ubuntu (latest), R release and devel

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
