// sblock.d
// Class for structured blocks of cells, for use within Eilmer4.
// This is the "classic" block within the mbcns/Eilmer series 
// of flow simulation codes.

// Peter J. 2014-07-20 first cut.

module sblock;

import std.conv;
import std.file;
import std.json;
import block;

class SBlock: Block {
public:
    int ni;
    int nj;
    int nk;

    this(int id, in char[] file_name) 
    {
	this.id = id;
	if ( file_name.length > 0 ) {
	    auto text = cast(string) read(file_name);
	    auto items = parseJSON(text);
	    ni = to!int(items["ni"].integer);
	    nj = to!int(items["nj"].integer);
	    nk = to!int(items["nk"].integer);
	}
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "SBlock(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ", ni=" ~ to!string(ni);
	repr ~= ", nj=" ~ to!string(nj);
	repr ~= ", nk=" ~ to!string(nk);
	repr ~= ")";
	return to!string(repr);
    }

private:

} // end class SBlock
