#! /bin/bash
# post-surface.sh
e3post.py --job=convex-ramp --tindx=last \
    --heat-flux-list="0,2,:,0,0;2,2,:,0,0;4,2,:,0,0;6,2,:,0,0;8,2,:,0,0;10,2,:,0,0;12,2,:,0,0;14,2,:,0,0;16,2,:,0,0;18,2,:,0,0;20,2,:,0,0;22,2,:,0,0;24,2,:,0,0;26,2,:,0,0" \
    --output-file=my-surface.data

echo "At this point, we should have data to view"


