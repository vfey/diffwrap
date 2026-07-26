---
title: CRAN package diffwrap
---

## Resubmission 2026-02-25
This is a resubmission of the package after addressing the problems found at initial submission,
namely in the Date, Description and License fields in the DESCRIPTION file.
The version was increased to 0.5-12.

## New submission 2026-02-25
I am submitting the latest stable version of the package which was finalized on Nov 7 2025.
Submission was delayed due to necessary modifications in strong dependency 'convertid'.

### Notes
The win-builder check returned a note about many dependencies. Given the purpose of 'diffwrap' as a wrapper package
for differential expression analysis a large number of dependencies and imported packages is unavoidable.

## Test environments (2026-02-25 - )
* local OS X install: x86_64-apple-darwin25.3.0, R 4.5.2
* win-builder (devel, release and oldrelease)
* Red Hat Enterprise Linux release 9.7 (Plow), R 4.5.2
