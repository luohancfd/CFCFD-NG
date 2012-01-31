#! /bin/bash
DEST=peterj@viking.mech.uq.edu.au:/var/www/html/cfcfd-new/

echo "About to rsync html files to viking."
echo "A password may have to be supplied for each rsync."
echo "Note that we assume several external directories"
echo "(of theses, etc) are available."

echo "1. Sphinx-generated HTML docs..."
rsync -av _build/html/ ${DEST}
echo "2. Theses and reports..."
rsync -av ../pdf ${DEST}
echo "3. Theses..."
rsync -av ~/papers/theses ${DEST}
echo "4. Case study GIFs..."
rsync -av ~/papers/drummond_simulation ${DEST}
echo "5. IMOC HTML docs..."
rsync -av ~/cfcfd3/app/imoc/doc/ ${DEST}/imoc/
echo "Done."


