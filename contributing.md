Contributing to this project:

- Any and all help is very much appreciated!
- Please make sure your branch is based off `stable`;
- Each PR should ideally add or change one feature only, or at most, a few closely related features;
- All PRs that change the C++ classes or the primarydock app *must* pass the master unit test (`tests/unit_test_master.sh`) before merge.
- All debug #defines must be set to zero (there are a great many of these at the end of `src/classes/constants.h`).
- The CFLAGS setting in the makefile must be set to release mode. We recommend setting it to debug before doing development, and then back to release before submitting the PR for review.
- Use whichever { style you prefer; as long as the code is readable and it works, that's all that should matter anyway.
- Have fun and try not to let the project vex you. (:
