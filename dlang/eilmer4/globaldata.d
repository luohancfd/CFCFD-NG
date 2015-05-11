/**
 * globaldata.d
 *
 * Author: Peter J. and Rowan G.
 * Versions: 2014-07-18 : first cut.
 *           2015-22-04 : added containers for solid blocks
 */

module globaldata;

import globalconfig;
import sblock;
import ssolidblock;

static SBlock[] allBlocks; // The array of Block objects, holding arrays of cells.
static SBlock[] myBlocks;  // Local collection that we can iterate over in parallel.

static SSolidBlock[] allSolidBlocks;
static SSolidBlock[] mySolidBlocks;

// The current (shared-memory) parallel code is based on having one SBlock per thread.
// We need to hava a dedicated set of configuration parameters for each thread so that
// there is no need to have memory barriers guarding their access.
static LocalConfig[] dedicatedConfig;
