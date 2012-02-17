#! /bin/bash

echo "About to rsync html files to viking."

echo "Please enter your username on viking:"
read username

echo "A password may have to be supplied for each rsync."
echo "Note that we assume several external directories"
echo "(of theses, etc) are available."

DEST=${username}@viking.mech.uq.edu.au:/var/www/html/cfcfd/

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
echo "6. Doxygen-generated HTML docs..."
rsync -av ../doxygen/html/ ${DEST}/doxygen/
echo "Done."


