#!/bin/bash
cd "$(dirname "$0")"
cd ..

clear

if make primarydock; then
    php -f predict/method_combined.php lig=2-methylbutylamine
    php -f predict/method_combined.php lig=N-methylpiperidine
    php -f predict/method_combined.php lig=ambroxide
    php -f predict/method_combined.php lig=beta-phenethylamine
    php -f predict/method_combined.php lig=coumarin
    php -f predict/method_combined.php lig=cyclohexylamine
    php -f predict/method_combined.php lig=eugenol
    php -f predict/method_combined.php lig=isoamylamine
    php -f predict/method_combined.php lig=timberol
    php -f predict/method_combined.php lig=triethylamine
    php -f predict/method_combined.php lig=trimethylamine
    php -f predict/method_combined.php lig=tyramine

    php -f predict/method_combined.php prot=TAAR6 lig=putrescine
    php -f predict/method_combined.php prot=TAAR6 lig=cadaverine

    php -f predict/method_combined.php prot=TAAR8 lig=putrescine
    php -f predict/method_combined.php prot=TAAR8 lig=cadaverine
    php -f predict/method_combined.php prot=TAAR8 lig=N-methylpiperidine

    php -f predict/method_combined.php prot=TAAR9 lig=N-methylpiperidine
fi

