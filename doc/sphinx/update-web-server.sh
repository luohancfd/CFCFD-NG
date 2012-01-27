#! /bin/bash
echo "About to rsync html files to viking"
echo "A password will have to be supplied for each rsync."
echo "First one..."
rsync -av _build/html/ peterj@viking.mech.uq.edu.au:/var/www/html/cfcfd-new/
echo "Second one..."
rsync -av ../pdf peterj@viking.mech.uq.edu.au:/var/www/html/cfcfd-new/
echo "Done."


