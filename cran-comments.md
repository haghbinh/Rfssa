## Resubmission of an archived package

Rfssa 3.2.0 was archived on CRAN on 2026-04-28 "for policy violation.
On Internet access." This version (3.2.1) fixes that problem:

* `loadCallcenterData()`, `loadJambiData()`, `loadMontanaData()`,
  `loadAustinData()` and `loadUtqiagvikData()` download raw example data from
  GitHub. They now fail gracefully, as the CRAN policy requires. If the resource
  is unavailable (no Internet connection, timeout, or a changed resource), they
  return `NULL` invisibly with an informative message instead of an error or
  warning.
* All examples that use these functions are wrapped in `\donttest{}` and check
  for a `NULL` result. The `funts()` example also has an offline example with
  simulated data.
* The Shiny apps started by `launchApp()` show a message instead of stopping
  when the data cannot be downloaded.
* `plot.funts()` now restores the user's `par()` settings on exit.

We checked the fix by running `R CMD check --as-cran` with the network blocked
(all `http(s)_proxy` variables pointing to a closed local port). Version 3.2.0
fails in the `funts` example ("cannot open URL ... Callcenter.rds"); version
3.2.1 passes, including the `--run-donttest` examples.

## Test environments

* Rocky Linux 9.8, R 4.6.1, GCC 11.5.0, reference BLAS (FlexiBLAS NETLIB),
  `OMP_THREAD_LIMIT=2`: `R CMD check --as-cran`.
* Same machine with its default multi-threaded OpenBLAS: `R CMD check --as-cran`
  with and without Internet access.

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking CRAN incoming feasibility ... NOTE
  New submission. Package was archived on CRAN (X-CRAN-Comment: Archived on
  2026-04-28 for policy violation. On Internet access.)
  This resubmission fixes the reason for the archival (see above).

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
  '-Werror=format-security' '-Wp,-D_FORTIFY_SOURCE=2'
  '-Wp,-D_GLIBCXX_ASSERTIONS' '-march=x86-64-v2'
  These flags come from the system R configuration (Makeconf) of the test
  machine, not from the package's `src/Makevars`.

* checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s:
  fssa (5.7 s user, 5.9 s elapsed). This example is unchanged from 3.2.0.

## Reverse dependencies

None known (the package was archived).
