Contributing to this project:

- I am sooo very grateful for the help!
- Please make sure your branch is based off `stable`;
- Each PR should ideally add or change one feature only, or at most, a few closely related features;
- All PRs that change the C++ classes or the primarydock app *must* pass the Big Three test (see the test/big_three file), and receive a three green squares message, before merge.
- All debug #defines must be set to zero (there are a great many of these at the end of `src/classes/constants.h`).
- The CFLAGS setting in the makefile must be set to release mode. We recommend setting it to debug before doing development, and then back to release before submitting the PR for review.
- Use whichever { style you prefer; as long as the code is readable and it works, that's all that should matter anyway.
- Have fun and try not to let the project vex you. (:
